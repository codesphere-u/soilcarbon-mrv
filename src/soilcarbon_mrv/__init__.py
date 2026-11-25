"""Equivalent Soil Mass (ESM) soil organic carbon MRV & uncertainty engine."""

from .certificate import MrvCertificate, issue_certificate, verify_certificate
from .esm import (
    CO2E_PER_C,
    EsmComparison,
    HorizonSample,
    InsufficientSoilMassError,
    SoilCore,
    compare_fixed_depth_vs_esm,
    compute_esm_stock,
    compute_fixed_depth_stock,
    pchip_evaluate,
    reference_soil_mass,
)
from .mrv import (
    PairedPlot,
    StratumMrvResult,
    evaluate_stratum_mrv,
    student_t_critical_90,
)

__all__ = [
    "CO2E_PER_C",
    "EsmComparison",
    "HorizonSample",
    "InsufficientSoilMassError",
    "MrvCertificate",
    "PairedPlot",
    "SoilCore",
    "StratumMrvResult",
    "compare_fixed_depth_vs_esm",
    "compute_esm_stock",
    "compute_fixed_depth_stock",
    "evaluate_stratum_mrv",
    "issue_certificate",
    "pchip_evaluate",
    "reference_soil_mass",
    "student_t_critical_90",
    "verify_certificate",
]
__version__ = "0.1.0"
