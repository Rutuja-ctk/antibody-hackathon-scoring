from dataclasses import dataclass
from typing import Optional, Dict, Any


@dataclass
class DesignInput:
    team_name: str
    design_id: str


@dataclass
class RawMetrics:
    # Binding and interface
    delta_g: Optional[float] = None
    contacts: Optional[int] = None
    dockq: Optional[float] = None
    iface_plddt: Optional[float] = None
    ipsae: Optional[float] = None
    cdr_sasa: Optional[float] = None

    # Developability
    netsolp: Optional[float] = None

    # Novelty
    cdr3_identity: Optional[float] = None
    anarci_pass: Optional[bool] = None


@dataclass
class PerMetricScores:
    # per metric scores on 0 to 10
    delta_g_score: float = 0.0
    contacts_score: float = 0.0
    dockq_score: float = 0.0
    ipsae_score: float = 0.0
    iface_plddt_score: float = 0.0
    cdr_sasa_score: float = 0.0
    netsolp_score: float = 0.0
    novelty_score: float = 0.0


@dataclass
class CategoryScores:
    binding_struct_score: float = 0.0
    developability_score: float = 0.0
    novelty_category_score: float = 0.0


@dataclass
class FinalScores:
    per_metric: PerMetricScores
    category: CategoryScores
    final_score_10: float
    final_score_100: float
    is_viable: bool
    fail_reason: Optional[str] = None


def row_from_raw_and_scores(
    design: DesignInput,
    raw: RawMetrics,
    scores: FinalScores,
) -> Dict[str, Any]:
    """
    Build a flat dict suitable for a pandas DataFrame row.
    """
    row = {
        "team_name": design.team_name,
        "design_id": design.design_id,
        # raw metrics
        "delta_g": raw.delta_g,
        "contacts": raw.contacts,
        "dockq": raw.dockq,
        "iface_plddt": raw.iface_plddt,
        "ipsae": raw.ipsae,
        "cdr_sasa": raw.cdr_sasa,
        "netsolp": raw.netsolp,
        "cdr3_identity": raw.cdr3_identity,
        "anarci_pass": raw.anarci_pass,
        # per metric scores
        "delta_g_score": scores.per_metric.delta_g_score,
        "contacts_score": scores.per_metric.contacts_score,
        "dockq_score": scores.per_metric.dockq_score,
        "ipsae_score": scores.per_metric.ipsae_score,
        "iface_plddt_score": scores.per_metric.iface_plddt_score,
        "cdr_sasa_score": scores.per_metric.cdr_sasa_score,
        "netsolp_score": scores.per_metric.netsolp_score,
        "novelty_score": scores.per_metric.novelty_score,
        # category scores
        "binding_struct_score": scores.category.binding_struct_score,
        "developability_score": scores.category.developability_score,
        "novelty_category_score": scores.category.novelty_category_score,
        # final
        "final_score_10": scores.final_score_10,
        "final_score_100": scores.final_score_100,
        "is_viable": scores.is_viable,
        "fail_reason": scores.fail_reason,
    }
    return row
