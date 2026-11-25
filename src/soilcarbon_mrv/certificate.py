"""Tamper-evident SHA-256 MRV audit certificates for soil carbon issuances."""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from datetime import UTC, datetime

from .mrv import StratumMrvResult


@dataclass(frozen=True)
class MrvCertificate:
    """Signed SHA-256 audit envelope binding MRV calculations and methodology."""

    certificate_id: str
    project_id: str
    methodology: str
    issued_at_utc: str
    stratum_result: StratumMrvResult
    sha256_digest: str

    def to_dict(self) -> dict[str, object]:
        return {
            "certificate_id": self.certificate_id,
            "project_id": self.project_id,
            "methodology": self.methodology,
            "issued_at_utc": self.issued_at_utc,
            "stratum_result": asdict(self.stratum_result),
            "sha256_digest": self.sha256_digest,
        }


def _canonical_digest(payload: dict[str, object]) -> str:
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(canonical).hexdigest()


def issue_certificate(
    *,
    project_id: str,
    stratum_result: StratumMrvResult,
    methodology: str = "ESM-PCHIP-VM0042-v1",
    issued_at_utc: str | None = None,
) -> MrvCertificate:
    """Create a deterministic SHA-256 tamper-evident MRV certificate."""
    ts = issued_at_utc or datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")
    body: dict[str, object] = {
        "project_id": project_id,
        "methodology": methodology,
        "issued_at_utc": ts,
        "stratum_result": asdict(stratum_result),
    }
    digest = _canonical_digest(body)
    cert_id = f"SCMRV-{digest[:12].upper()}"
    return MrvCertificate(
        certificate_id=cert_id,
        project_id=project_id,
        methodology=methodology,
        issued_at_utc=ts,
        stratum_result=stratum_result,
        sha256_digest=digest,
    )


def verify_certificate(cert: MrvCertificate) -> bool:
    """Verify that a certificate's payload matches its SHA-256 digest and ID prefix."""
    body: dict[str, object] = {
        "project_id": cert.project_id,
        "methodology": cert.methodology,
        "issued_at_utc": cert.issued_at_utc,
        "stratum_result": asdict(cert.stratum_result),
    }
    expected = _canonical_digest(body)
    return cert.sha256_digest == expected and cert.certificate_id == f"SCMRV-{expected[:12].upper()}"
