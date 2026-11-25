"""Paired-plot stratum SOC accounting with Verra VM0042-style uncertainty & buffer deductions."""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass

from .esm import CO2E_PER_C, SoilCore, compute_esm_stock, compute_fixed_depth_stock, reference_soil_mass

# Two-tailed 90% CI (95% one-sided, alpha=0.05) Student's t critical values for df=1..30
_T_CRIT_90: dict[int, float] = {
    1: 6.3138,
    2: 2.9200,
    3: 2.3534,
    4: 2.1318,
    5: 2.0150,
    6: 1.9432,
    7: 1.8946,
    8: 1.8595,
    9: 1.8331,
    10: 1.8125,
    11: 1.7959,
    12: 1.7823,
    13: 1.7709,
    14: 1.7613,
    15: 1.7531,
    16: 1.7459,
    17: 1.7396,
    18: 1.7341,
    19: 1.7291,
    20: 1.7247,
    21: 1.7207,
    22: 1.7171,
    23: 1.7139,
    24: 1.7109,
    25: 1.7081,
    26: 1.7056,
    27: 1.7033,
    28: 1.7011,
    29: 1.6991,
    30: 1.6973,
}


def student_t_critical_90(df: int) -> float:
    """Return Student's t critical value for 90% two-tailed (95% one-sided) CI."""
    if df < 1:
        raise ValueError("Degrees of freedom must be >= 1")
    if df in _T_CRIT_90:
        return _T_CRIT_90[df]
    # Cornish-Fisher expansion for z = 1.6448536269514722 (accurate to <1e-4 for df > 30)
    z = 1.6448536269514722
    z3 = z**3
    z5 = z**5
    return z + (z3 + z) / (4.0 * df) + (5.0 * z5 + 16.0 * z3 + 3.0 * z) / (96.0 * df * df)


@dataclass(frozen=True)
class PairedPlot:
    """Baseline and monitoring soil cores sampled at the same plot."""

    plot_id: str
    baseline: SoilCore
    monitoring: SoilCore


@dataclass(frozen=True)
class StratumMrvResult:
    """Complete MRV carbon accounting result for a project stratum."""

    stratum_id: str
    area_ha: float
    plot_count: int
    reference_mass_mg_ha: float
    mean_delta_esm_mg_c_ha: float
    mean_delta_fd_mg_c_ha: float
    phantom_fd_bias_tco2e_total: float
    std_delta_tco2e_ha: float
    se_delta_tco2e_ha: float
    t_critical_90: float
    relative_uncertainty_pct: float
    gross_removals_tco2e: float
    uncertainty_deduction_pct: float
    uncertainty_deduction_tco2e: float
    buffer_pool_pct: float
    buffer_pool_tco2e: float
    leakage_pct: float
    leakage_tco2e: float
    net_issuable_tco2e: float
    statistically_significant: bool


def evaluate_stratum_mrv(
    stratum_id: str,
    area_ha: float,
    plots: Sequence[PairedPlot],
    *,
    target_depth_cm: float = 30.0,
    reference_mass_mg_ha: float | None = None,
    uncertainty_threshold_pct: float = 15.0,
    buffer_pool_pct: float = 20.0,
    leakage_pct: float = 0.0,
    allow_linear_extrapolation: bool = False,
) -> StratumMrvResult:
    """Evaluate paired-plot ESM carbon stock change and conservative credit deductions."""
    if area_ha <= 0.0:
        raise ValueError("area_ha must be positive")
    if len(plots) < 2:
        raise ValueError("At least 2 paired plots are required for stratum uncertainty analysis")
    if not (0.0 <= buffer_pool_pct < 100.0):
        raise ValueError("buffer_pool_pct must be in [0, 100)")
    if not (0.0 <= leakage_pct < 100.0):
        raise ValueError("leakage_pct must be in [0, 100)")

    all_baseline = [p.baseline for p in plots]
    ref_mass = (
        reference_mass_mg_ha
        if reference_mass_mg_ha is not None
        else reference_soil_mass(all_baseline, target_depth_cm)
    )

    deltas_esm_c: list[float] = []
    deltas_fd_c: list[float] = []
    for plot in plots:
        b_esm = compute_esm_stock(
            plot.baseline,
            ref_mass,
            allow_linear_extrapolation=allow_linear_extrapolation,
        )
        m_esm = compute_esm_stock(
            plot.monitoring,
            ref_mass,
            allow_linear_extrapolation=allow_linear_extrapolation,
        )
        deltas_esm_c.append(m_esm - b_esm)

        b_fd = compute_fixed_depth_stock(plot.baseline, target_depth_cm)
        m_fd = compute_fixed_depth_stock(plot.monitoring, target_depth_cm)
        deltas_fd_c.append(m_fd - b_fd)

    n = len(plots)
    mean_esm_c = sum(deltas_esm_c) / n
    mean_fd_c = sum(deltas_fd_c) / n

    deltas_co2e = [d * CO2E_PER_C for d in deltas_esm_c]
    mean_co2e_ha = mean_esm_c * CO2E_PER_C
    variance = sum((x - mean_co2e_ha) ** 2 for x in deltas_co2e) / (n - 1)
    std_co2e_ha = math.sqrt(variance)
    se_co2e_ha = std_co2e_ha / math.sqrt(n)
    t_crit = student_t_critical_90(n - 1)
    half_width = t_crit * se_co2e_ha

    gross_tco2e = mean_co2e_ha * area_ha
    phantom_bias_total = (mean_fd_c - mean_esm_c) * CO2E_PER_C * area_ha

    if mean_co2e_ha <= 0.0:
        rel_unc_pct = 100.0
        unc_deduct_pct = 100.0
        unc_deduct_tco2e = 0.0
        buffer_tco2e = 0.0
        leakage_tco2e = 0.0
        net_issuable = 0.0
        significant = False
    else:
        rel_unc_pct = (half_width / mean_co2e_ha) * 100.0
        significant = (mean_co2e_ha - half_width) > 0.0
        # VM0042-style deduction: deduct excess uncertainty above threshold (default 15%)
        unc_deduct_pct = min(100.0, max(0.0, rel_unc_pct - uncertainty_threshold_pct))
        unc_deduct_tco2e = gross_tco2e * (unc_deduct_pct / 100.0)
        after_unc = max(0.0, gross_tco2e - unc_deduct_tco2e)
        leakage_tco2e = after_unc * (leakage_pct / 100.0)
        after_leakage = max(0.0, after_unc - leakage_tco2e)
        buffer_tco2e = after_leakage * (buffer_pool_pct / 100.0)
        net_issuable = max(0.0, after_leakage - buffer_tco2e)

    return StratumMrvResult(
        stratum_id=stratum_id,
        area_ha=round(area_ha, 4),
        plot_count=n,
        reference_mass_mg_ha=round(ref_mass, 2),
        mean_delta_esm_mg_c_ha=round(mean_esm_c, 4),
        mean_delta_fd_mg_c_ha=round(mean_fd_c, 4),
        phantom_fd_bias_tco2e_total=round(phantom_bias_total, 2),
        std_delta_tco2e_ha=round(std_co2e_ha, 4),
        se_delta_tco2e_ha=round(se_co2e_ha, 4),
        t_critical_90=round(t_crit, 4),
        relative_uncertainty_pct=round(rel_unc_pct, 2),
        gross_removals_tco2e=round(gross_tco2e, 2),
        uncertainty_deduction_pct=round(unc_deduct_pct, 2),
        uncertainty_deduction_tco2e=round(unc_deduct_tco2e, 2),
        buffer_pool_pct=round(buffer_pool_pct, 2),
        buffer_pool_tco2e=round(buffer_tco2e, 2),
        leakage_pct=round(leakage_pct, 2),
        leakage_tco2e=round(leakage_tco2e, 2),
        net_issuable_tco2e=round(net_issuable, 2),
        statistically_significant=significant,
    )
