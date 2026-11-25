"""Tests for soilcarbon-mrv."""

from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path

import pytest

from soilcarbon_mrv import (
    HorizonSample,
    InsufficientSoilMassError,
    PairedPlot,
    SoilCore,
    compare_fixed_depth_vs_esm,
    compute_esm_stock,
    compute_fixed_depth_stock,
    evaluate_stratum_mrv,
    issue_certificate,
    pchip_evaluate,
    reference_soil_mass,
    student_t_critical_90,
    verify_certificate,
)
from soilcarbon_mrv.cli import main as cli_main


def _make_core(core_id: str, bd_top: float, soc_top: float, bd_sub: float, soc_sub: float) -> SoilCore:
    return SoilCore(
        core_id=core_id,
        horizons=(
            HorizonSample(0.0, 15.0, bd_top, soc_top),
            HorizonSample(15.0, 30.0, bd_sub, soc_sub),
            HorizonSample(30.0, 45.0, bd_sub, soc_sub * 0.7),
        ),
    )


def test_compaction_phantom_carbon_eliminated_by_esm() -> None:
    """When bulk density rises with identical SOC%, fixed depth hallucinates gains while ESM is ~0."""
    baseline = SoilCore(
        core_id="B1",
        horizons=(
            HorizonSample(0.0, 15.0, 1.20, 2.0),
            HorizonSample(15.0, 30.0, 1.20, 2.0),
        ),
    )
    compacted = SoilCore(
        core_id="M1",
        horizons=(
            HorizonSample(0.0, 15.0, 1.35, 2.0),
            HorizonSample(15.0, 30.0, 1.35, 2.0),
        ),
    )
    cmp = compare_fixed_depth_vs_esm(baseline, compacted, target_depth_cm=30.0)
    # Fixed depth falsely claims +9.0 Mg C/ha (+33.0 tCO2e/ha)
    assert cmp.delta_fd_mg_c_ha == pytest.approx(9.0, abs=1e-3)
    # ESM correctly sees zero true carbon increase at the 3600 Mg/ha reference mass
    assert cmp.delta_esm_mg_c_ha == pytest.approx(0.0, abs=1e-4)
    assert cmp.phantom_fd_bias_tco2e_ha == pytest.approx(33.0, abs=0.05)


def test_coarse_fragment_reduces_fine_earth_mass() -> None:
    rocky = HorizonSample(0.0, 30.0, 1.30, 1.5, coarse_fragment_pct=20.0)
    assert rocky.fine_earth_mass_mg_ha == pytest.approx(1.30 * 30.0 * 100.0 * 0.80)
    assert rocky.soc_stock_mg_c_ha == pytest.approx(rocky.fine_earth_mass_mg_ha * 0.015)


def test_insufficient_soil_mass_raises_unless_extrapolate_enabled() -> None:
    shallow = SoilCore("S1", (HorizonSample(0.0, 20.0, 1.1, 1.8),))
    with pytest.raises(InsufficientSoilMassError):
        compute_esm_stock(shallow, reference_mass_mg_ha=3000.0)
    extrap = compute_esm_stock(shallow, reference_mass_mg_ha=3000.0, allow_linear_extrapolation=True)
    assert extrap == pytest.approx(3000.0 * 0.018)


def test_pchip_monotonicity_preserved() -> None:
    xs = [0.0, 1500.0, 3500.0, 5500.0]
    ys = [0.0, 42.0, 62.0, 70.0]
    prev = -1.0
    for q in range(0, 5501, 250):
        val = pchip_evaluate(xs, ys, float(q))
        assert val >= prev - 1e-9
        prev = val


def test_stratum_mrv_deductions_and_certificate_integrity() -> None:
    plots = [
        PairedPlot("P1", _make_core("B1", 1.20, 1.60, 1.30, 0.90), _make_core("M1", 1.24, 1.92, 1.32, 0.98)),
        PairedPlot("P2", _make_core("B2", 1.18, 1.55, 1.28, 0.85), _make_core("M2", 1.22, 1.88, 1.30, 0.93)),
        PairedPlot("P3", _make_core("B3", 1.22, 1.65, 1.31, 0.92), _make_core("M3", 1.25, 1.96, 1.33, 1.00)),
        PairedPlot("P4", _make_core("B4", 1.19, 1.58, 1.29, 0.88), _make_core("M4", 1.23, 1.90, 1.31, 0.95)),
    ]
    res = evaluate_stratum_mrv("EAST-UG-S1", area_ha=250.0, plots=plots, buffer_pool_pct=20.0, leakage_pct=5.0)
    assert res.statistically_significant is True
    assert res.gross_removals_tco2e > 0.0
    assert res.relative_uncertainty_pct < 15.0
    assert res.uncertainty_deduction_pct == 0.0
    assert res.net_issuable_tco2e < res.gross_removals_tco2e

    cert = issue_certificate(project_id="UG-REGEN-2026", stratum_result=res, issued_at_utc="2026-11-25T12:00:00Z")
    assert verify_certificate(cert) is True

    tampered = replace(cert, stratum_result=replace(res, net_issuable_tco2e=res.net_issuable_tco2e + 100.0))
    assert verify_certificate(tampered) is False


def test_high_uncertainty_triggers_vm0042_deduction() -> None:
    plots = [
        PairedPlot("P1", _make_core("B1", 1.20, 1.50, 1.30, 0.90), _make_core("M1", 1.20, 2.20, 1.30, 0.95)),
        PairedPlot("P2", _make_core("B2", 1.20, 1.50, 1.30, 0.90), _make_core("M2", 1.20, 1.55, 1.30, 0.91)),
        PairedPlot("P3", _make_core("B3", 1.20, 1.50, 1.30, 0.90), _make_core("M3", 1.20, 1.60, 1.30, 0.92)),
    ]
    res = evaluate_stratum_mrv("NOISY-STRATUM", area_ha=100.0, plots=plots)
    assert res.relative_uncertainty_pct > 15.0
    assert res.uncertainty_deduction_pct > 0.0
    assert res.uncertainty_deduction_tco2e > 0.0


def test_student_t_high_df_and_reference_mass() -> None:
    assert student_t_critical_90(100) == pytest.approx(1.6602, abs=0.002)
    c1 = _make_core("C1", 1.15, 1.5, 1.25, 0.8)
    c2 = _make_core("C2", 1.25, 1.5, 1.35, 0.8)
    assert reference_soil_mass([c1, c2], 30.0) == pytest.approx(3600.0)
    assert compute_fixed_depth_stock(c1, 15.0) == pytest.approx(1.15 * 15 * 100 * 0.015)


def test_cli_outputs_valid_certificate(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    payload = {
        "project_id": "DEMO-SOC",
        "stratum_id": "S1",
        "area_ha": 50.0,
        "plots": [
            {
                "plot_id": "P1",
                "baseline": {
                    "core_id": "B1",
                    "horizons": [
                        {"depth_top_cm": 0, "depth_bottom_cm": 15, "bulk_density_g_cm3": 1.2, "soc_percent": 1.5},
                        {"depth_top_cm": 15, "depth_bottom_cm": 35, "bulk_density_g_cm3": 1.3, "soc_percent": 0.9},
                    ],
                },
                "monitoring": {
                    "core_id": "M1",
                    "horizons": [
                        {"depth_top_cm": 0, "depth_bottom_cm": 15, "bulk_density_g_cm3": 1.22, "soc_percent": 1.85},
                        {"depth_top_cm": 15, "depth_bottom_cm": 35, "bulk_density_g_cm3": 1.32, "soc_percent": 0.95},
                    ],
                },
            },
            {
                "plot_id": "P2",
                "baseline": {
                    "core_id": "B2",
                    "horizons": [
                        {"depth_top_cm": 0, "depth_bottom_cm": 15, "bulk_density_g_cm3": 1.18, "soc_percent": 1.48},
                        {"depth_top_cm": 15, "depth_bottom_cm": 35, "bulk_density_g_cm3": 1.28, "soc_percent": 0.88},
                    ],
                },
                "monitoring": {
                    "core_id": "M2",
                    "horizons": [
                        {"depth_top_cm": 0, "depth_bottom_cm": 15, "bulk_density_g_cm3": 1.21, "soc_percent": 1.80},
                        {"depth_top_cm": 15, "depth_bottom_cm": 35, "bulk_density_g_cm3": 1.30, "soc_percent": 0.93},
                    ],
                },
            },
        ],
    }
    in_file = tmp_path / "stratum.json"
    in_file.write_text(json.dumps(payload), encoding="utf-8")
    assert cli_main([str(in_file)]) == 0
    out = json.loads(capsys.readouterr().out)
    assert out["certificate_id"].startswith("SCMRV-")
    assert out["stratum_result"]["net_issuable_tco2e"] > 0.0
