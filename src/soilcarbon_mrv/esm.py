"""Equivalent Soil Mass (ESM) soil organic carbon stock calculation.

Implements rock-fragment-corrected soil mass, monotone cubic Hermite (PCHIP)
interpolation of cumulative SOC vs. cumulative fine-earth mineral+organic mass
(Wendt & Hauser 2013; von Haden et al. 2020), and fixed-depth bias diagnostics.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

CO2E_PER_C: float = 44.0 / 12.0


class InsufficientSoilMassError(ValueError):
    """Raised when a core's sampled cumulative soil mass is below the reference mass."""


@dataclass(frozen=True)
class HorizonSample:
    """Single soil depth increment measured from a core."""

    depth_top_cm: float
    depth_bottom_cm: float
    bulk_density_g_cm3: float
    soc_percent: float
    coarse_fragment_pct: float = 0.0

    def __post_init__(self) -> None:
        if self.depth_top_cm < 0.0 or self.depth_bottom_cm <= self.depth_top_cm:
            raise ValueError(f"Invalid depth interval: {self.depth_top_cm}..{self.depth_bottom_cm} cm")
        if self.bulk_density_g_cm3 <= 0.0 or self.bulk_density_g_cm3 > 2.65:
            raise ValueError(f"Bulk density must be in (0, 2.65] g/cm3, got {self.bulk_density_g_cm3}")
        if self.soc_percent < 0.0 or self.soc_percent > 58.0:
            raise ValueError(f"SOC percent must be in [0, 58], got {self.soc_percent}")
        if self.coarse_fragment_pct < 0.0 or self.coarse_fragment_pct >= 100.0:
            raise ValueError(f"Coarse fragment percent must be in [0, 100), got {self.coarse_fragment_pct}")

    @property
    def thickness_cm(self) -> float:
        return self.depth_bottom_cm - self.depth_top_cm

    @property
    def fine_earth_mass_mg_ha(self) -> float:
        """Fine-earth (<2mm) soil mass in Mg/ha (t/ha).

        1 g/cm3 * 1 cm = 1 g/cm2 = 100 Mg/ha.
        """
        fine_fraction = 1.0 - (self.coarse_fragment_pct / 100.0)
        return self.bulk_density_g_cm3 * self.thickness_cm * 100.0 * fine_fraction

    @property
    def soc_stock_mg_c_ha(self) -> float:
        """Horizon organic carbon stock in Mg C/ha."""
        return self.fine_earth_mass_mg_ha * (self.soc_percent / 100.0)


@dataclass(frozen=True)
class SoilCore:
    """Ordered sequence of contiguous depth increments from a single sampling point."""

    core_id: str
    horizons: tuple[HorizonSample, ...]
    latitude: float | None = None
    longitude: float | None = None

    def __post_init__(self) -> None:
        if not self.core_id.strip():
            raise ValueError("core_id must be non-empty")
        if not self.horizons:
            raise ValueError("SoilCore must contain at least one HorizonSample")
        if abs(self.horizons[0].depth_top_cm) > 1e-6:
            raise ValueError(f"First horizon must start at 0 cm, got {self.horizons[0].depth_top_cm}")
        for prev, curr in zip(self.horizons, self.horizons[1:], strict=False):
            if abs(curr.depth_top_cm - prev.depth_bottom_cm) > 1e-6:
                raise ValueError(
                    f"Non-contiguous horizons in core {self.core_id}: "
                    f"{prev.depth_bottom_cm} cm -> {curr.depth_top_cm} cm"
                )

    @property
    def max_depth_cm(self) -> float:
        return self.horizons[-1].depth_bottom_cm

    @property
    def total_mass_mg_ha(self) -> float:
        return sum(h.fine_earth_mass_mg_ha for h in self.horizons)

    @property
    def total_soc_mg_c_ha(self) -> float:
        return sum(h.soc_stock_mg_c_ha for h in self.horizons)

    def cumulative_profile(self) -> tuple[list[float], list[float]]:
        """Return cumulative (fine_earth_mass_mg_ha, soc_stock_mg_c_ha) starting from (0, 0)."""
        masses = [0.0]
        carbons = [0.0]
        cum_m = 0.0
        cum_c = 0.0
        for h in self.horizons:
            cum_m += h.fine_earth_mass_mg_ha
            cum_c += h.soc_stock_mg_c_ha
            masses.append(cum_m)
            carbons.append(cum_c)
        return masses, carbons


def pchip_evaluate(xs: Sequence[float], ys: Sequence[float], x_query: float) -> float:
    """Evaluate monotone Piecewise Cubic Hermite Interpolating Polynomial (Fritsch-Carlson).

    Preserves monotonicity so cumulative SOC stock cannot overshoot or decrease.
    """
    n = len(xs)
    if n != len(ys) or n < 2:
        raise ValueError("pchip_evaluate requires matching sequences of length >= 2")
    if x_query < xs[0] or x_query > xs[-1]:
        raise ValueError(f"x_query={x_query} out of range [{xs[0]}, {xs[-1]}]")
    if x_query == xs[0]:
        return ys[0]
    if x_query == xs[-1]:
        return ys[-1]

    hs = [xs[i + 1] - xs[i] for i in range(n - 1)]
    if any(h <= 0.0 for h in hs):
        raise ValueError("xs must be strictly increasing")
    deltas = [(ys[i + 1] - ys[i]) / hs[i] for i in range(n - 1)]

    if n == 2:
        return ys[0] + deltas[0] * (x_query - xs[0])

    # Compute Fritsch-Carlson slopes d[i]
    d = [0.0] * n
    for i in range(1, n - 1):
        if deltas[i - 1] == 0.0 or deltas[i] == 0.0 or (deltas[i - 1] * deltas[i] < 0.0):
            d[i] = 0.0
        else:
            w1 = 2.0 * hs[i] + hs[i - 1]
            w2 = hs[i] + 2.0 * hs[i - 1]
            d[i] = (w1 + w2) / (w1 / deltas[i - 1] + w2 / deltas[i])

    # Endpoint slopes (shape-preserving three-point formula)
    d0 = ((2.0 * hs[0] + hs[1]) * deltas[0] - hs[0] * deltas[1]) / (hs[0] + hs[1])
    if d0 * deltas[0] <= 0.0:
        d0 = 0.0
    elif deltas[0] * deltas[1] < 0.0 and abs(d0) > 3.0 * abs(deltas[0]):
        d0 = 3.0 * deltas[0]
    d[0] = d0

    dn = ((2.0 * hs[-1] + hs[-2]) * deltas[-1] - hs[-1] * deltas[-2]) / (hs[-1] + hs[-2])
    if dn * deltas[-1] <= 0.0:
        dn = 0.0
    elif deltas[-1] * deltas[-2] < 0.0 and abs(dn) > 3.0 * abs(deltas[-1]):
        dn = 3.0 * deltas[-1]
    d[-1] = dn

    # Locate interval
    idx = 0
    for i in range(n - 1):
        if xs[i] <= x_query <= xs[i + 1]:
            idx = i
            break

    h = hs[idx]
    t = (x_query - xs[idx]) / h
    t2 = t * t
    t3 = t2 * t
    h00 = 2.0 * t3 - 3.0 * t2 + 1.0
    h10 = t3 - 2.0 * t2 + t
    h01 = -2.0 * t3 + 3.0 * t2
    h11 = t3 - t2
    return h00 * ys[idx] + h10 * h * d[idx] + h01 * ys[idx + 1] + h11 * h * d[idx + 1]


def compute_fixed_depth_stock(core: SoilCore, target_depth_cm: float) -> float:
    """Compute naive fixed-depth SOC stock (Mg C/ha) down to target_depth_cm."""
    if target_depth_cm <= 0.0:
        raise ValueError("target_depth_cm must be positive")
    if target_depth_cm > core.max_depth_cm + 1e-6:
        raise ValueError(f"Core {core.core_id} only reaches {core.max_depth_cm} cm (< {target_depth_cm} cm)")

    total_c = 0.0
    for h in core.horizons:
        if h.depth_top_cm >= target_depth_cm:
            break
        used_thickness = min(h.depth_bottom_cm, target_depth_cm) - h.depth_top_cm
        frac = used_thickness / h.thickness_cm
        total_c += h.soc_stock_mg_c_ha * frac
    return total_c


def reference_soil_mass(cores: Sequence[SoilCore], target_depth_cm: float) -> float:
    """Select the minimum cumulative fine-earth mass (Mg/ha) at target_depth_cm across cores.

    Using the minimum baseline mass avoids extrapolation into unsampled subsoil.
    """
    if not cores:
        raise ValueError("At least one SoilCore is required to determine reference soil mass")
    masses: list[float] = []
    for core in cores:
        if target_depth_cm > core.max_depth_cm + 1e-6:
            raise ValueError(f"Core {core.core_id} depth {core.max_depth_cm} cm < target {target_depth_cm} cm")
        m = 0.0
        for h in core.horizons:
            if h.depth_top_cm >= target_depth_cm:
                break
            used = min(h.depth_bottom_cm, target_depth_cm) - h.depth_top_cm
            m += h.fine_earth_mass_mg_ha * (used / h.thickness_cm)
        masses.append(m)
    return min(masses)


def compute_esm_stock(
    core: SoilCore,
    reference_mass_mg_ha: float,
    *,
    allow_linear_extrapolation: bool = False,
) -> float:
    """Compute Equivalent Soil Mass SOC stock (Mg C/ha) at reference_mass_mg_ha."""
    if reference_mass_mg_ha <= 0.0:
        raise ValueError("reference_mass_mg_ha must be positive")
    masses, carbons = core.cumulative_profile()
    if reference_mass_mg_ha <= masses[-1] + 1e-6:
        clamped = min(reference_mass_mg_ha, masses[-1])
        return pchip_evaluate(masses, carbons, clamped)
    if not allow_linear_extrapolation:
        raise InsufficientSoilMassError(
            f"Core {core.core_id} total mass {masses[-1]:.1f} Mg/ha < reference {reference_mass_mg_ha:.1f} Mg/ha. "
            "Sample deeper or enable explicit extrapolation."
        )
    # Extrapolate using deepest horizon concentration
    excess_mass = reference_mass_mg_ha - masses[-1]
    deepest_soc_frac = core.horizons[-1].soc_percent / 100.0
    return carbons[-1] + excess_mass * deepest_soc_frac


@dataclass(frozen=True)
class EsmComparison:
    """Diagnostic comparing naive Fixed-Depth vs. Equivalent Soil Mass carbon delta."""

    baseline_fd_mg_c_ha: float
    monitoring_fd_mg_c_ha: float
    delta_fd_mg_c_ha: float
    baseline_esm_mg_c_ha: float
    monitoring_esm_mg_c_ha: float
    delta_esm_mg_c_ha: float
    phantom_fd_bias_mg_c_ha: float
    phantom_fd_bias_tco2e_ha: float
    reference_mass_mg_ha: float


def compare_fixed_depth_vs_esm(
    baseline: SoilCore,
    monitoring: SoilCore,
    target_depth_cm: float,
    reference_mass_mg_ha: float | None = None,
) -> EsmComparison:
    """Quantify how much phantom carbon gain/loss fixed-depth accounting introduces."""
    ref_mass = (
        reference_mass_mg_ha
        if reference_mass_mg_ha is not None
        else reference_soil_mass([baseline, monitoring], target_depth_cm)
    )
    base_fd = compute_fixed_depth_stock(baseline, target_depth_cm)
    mon_fd = compute_fixed_depth_stock(monitoring, target_depth_cm)
    delta_fd = mon_fd - base_fd

    base_esm = compute_esm_stock(baseline, ref_mass)
    mon_esm = compute_esm_stock(monitoring, ref_mass)
    delta_esm = mon_esm - base_esm

    bias_c = delta_fd - delta_esm
    return EsmComparison(
        baseline_fd_mg_c_ha=round(base_fd, 4),
        monitoring_fd_mg_c_ha=round(mon_fd, 4),
        delta_fd_mg_c_ha=round(delta_fd, 4),
        baseline_esm_mg_c_ha=round(base_esm, 4),
        monitoring_esm_mg_c_ha=round(mon_esm, 4),
        delta_esm_mg_c_ha=round(delta_esm, 4),
        phantom_fd_bias_mg_c_ha=round(bias_c, 4),
        phantom_fd_bias_tco2e_ha=round(bias_c * CO2E_PER_C, 4),
        reference_mass_mg_ha=round(ref_mass, 2),
    )
