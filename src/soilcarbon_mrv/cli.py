"""CLI entrypoint for soilcarbon-mrv."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from .certificate import issue_certificate
from .esm import HorizonSample, SoilCore
from .mrv import PairedPlot, evaluate_stratum_mrv


def _parse_core(raw: dict[str, Any]) -> SoilCore:
    horizons = tuple(
        HorizonSample(
            depth_top_cm=float(h["depth_top_cm"]),
            depth_bottom_cm=float(h["depth_bottom_cm"]),
            bulk_density_g_cm3=float(h["bulk_density_g_cm3"]),
            soc_percent=float(h["soc_percent"]),
            coarse_fragment_pct=float(h.get("coarse_fragment_pct", 0.0)),
        )
        for h in raw["horizons"]
    )
    return SoilCore(
        core_id=str(raw["core_id"]),
        horizons=horizons,
        latitude=float(raw["latitude"]) if raw.get("latitude") is not None else None,
        longitude=float(raw["longitude"]) if raw.get("longitude") is not None else None,
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="soilcarbon-mrv",
        description="Evaluate Equivalent Soil Mass (ESM) soil organic carbon MRV & uncertainty deductions.",
    )
    parser.add_argument("input_json", type=Path, help="Path to stratum paired-plot JSON payload.")
    parser.add_argument("--target-depth-cm", type=float, default=30.0, help="Target depth in cm (default: 30).")
    parser.add_argument("--buffer-pct", type=float, default=20.0, help="Permanence buffer pool percent (default: 20).")
    parser.add_argument("--leakage-pct", type=float, default=0.0, help="Leakage deduction percent (default: 0).")
    parser.add_argument("--extrapolate", action="store_true", help="Allow linear subsoil extrapolation if under-mass.")
    args = parser.parse_args(argv)

    data = json.loads(args.input_json.read_text(encoding="utf-8"))
    plots = [
        PairedPlot(
            plot_id=str(p["plot_id"]),
            baseline=_parse_core(p["baseline"]),
            monitoring=_parse_core(p["monitoring"]),
        )
        for p in data["plots"]
    ]
    result = evaluate_stratum_mrv(
        stratum_id=str(data["stratum_id"]),
        area_ha=float(data["area_ha"]),
        plots=plots,
        target_depth_cm=args.target_depth_cm,
        buffer_pool_pct=args.buffer_pct,
        leakage_pct=args.leakage_pct,
        allow_linear_extrapolation=args.extrapolate,
    )
    cert = issue_certificate(project_id=str(data.get("project_id", "PROJECT-01")), stratum_result=result)
    print(json.dumps(cert.to_dict(), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
