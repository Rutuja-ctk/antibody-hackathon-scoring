from typing import Optional, List
from .models import RawMetrics, PerMetricScores, CategoryScores, FinalScores


def _clamp(x: float, lo: float, hi: float) -> float:
    return max(lo, min(hi, x))


def _linear_scale(x: float, x0: float, x1: float, y0: float, y1: float) -> float:
    """
    Linearly map x in [x0, x1] to y in [y0, y1].
    Values outside [x0, x1] are clamped to y0 or y1.
    """
    if x0 == x1:
        return y0
    if x <= x0:
        return y0
    if x >= x1:
        return y1
    t = (x - x0) / (x1 - x0)
    return y0 + t * (y1 - y0)


# per metric scoring functions

def _score_ipsae(ipsae: Optional[float]) -> float:
    """
    ipSAE in [0,1], higher is better.

    Bands:
      0.00 - 0.40 -> 0
      0.40 - 0.60 -> 0 to 4   poor
      0.60 - 0.80 -> 4 to 8   medium
      0.80 - 1.00 -> 8 to 10  good
    """
    if ipsae is None:
        return 0.0
    x = _clamp(ipsae, 0.0, 1.0)

    if x <= 0.40:
        return 0.0
    elif x <= 0.60:
        return _linear_scale(x, 0.40, 0.60, 0.0, 4.0)
    elif x <= 0.80:
        return _linear_scale(x, 0.60, 0.80, 4.0, 8.0)
    else:
        return _linear_scale(x, 0.80, 1.00, 8.0, 10.0)


def _score_dockq(dockq: Optional[float]) -> float:
    """
    DockQ in [0,1], higher is better.

    CAPRI style bands:
      0.00 - 0.23   incorrect or poor
      0.23 - 0.49   acceptable
      0.49 - 0.80   medium
      0.80 - 1.00   high

    We map to:
      0.00 - 0.23  -> 0 to 3
      0.23 - 0.49  -> 3 to 5.5
      0.49 - 0.80  -> 5.5 to 8.5
      0.80 - 1.00  -> 8.5 to 10
    """
    if dockq is None:
        return 0.0
    x = _clamp(dockq, 0.0, 1.0)

    if x <= 0.23:
        return _linear_scale(x, 0.0, 0.23, 0.0, 3.0)
    elif x <= 0.49:
        return _linear_scale(x, 0.23, 0.49, 3.0, 5.5)
    elif x <= 0.80:
        return _linear_scale(x, 0.49, 0.80, 5.5, 8.5)
    else:
        return _linear_scale(x, 0.80, 1.00, 8.5, 10.0)


def _score_delta_g(delta_g: Optional[float]) -> float:
    """
    Delta G (kcal per mol), lower is better.

    Bands:
      >= -6        worst
      -6 to -10    poor
      -10 to -12   medium
      <= -12       good

    Map:
      >= -6         -> 0
      -6 to -10     -> 0 to 4
      -10 to -12    -> 4 to 8
      <= -12        -> 8 to 10
    """
    if delta_g is None:
        return 0.0
    x = delta_g

    if x >= -6.0:
        return 0.0
    elif x >= -10.0:
        return _linear_scale(x, -6.0, -10.0, 0.0, 4.0)
    elif x >= -12.0:
        return _linear_scale(x, -10.0, -12.0, 4.0, 8.0)
    else:
        return 10.0


def _score_contacts(contacts: Optional[int]) -> float:
    """
    Intermolecular contacts, higher is better.

    Bands:
      <= 10         worst
      10 - 15       poor
      15 - 25       medium
      >= 25         good
    """
    if contacts is None:
        return 0.0
    x = float(contacts)

    if x <= 10.0:
        return 0.0
    elif x <= 15.0:
        return _linear_scale(x, 10.0, 15.0, 0.0, 4.0)
    elif x <= 25.0:
        return _linear_scale(x, 15.0, 25.0, 4.0, 8.0)
    else:
        return _linear_scale(x, 25.0, 50.0, 8.0, 10.0)


def _score_iface_plddt(iface_plddt: Optional[float]) -> float:
    """
    Interface pLDDT, higher is better.

    Bands:
      < 65          poor
      65 - 70       poor to medium
      70 - 80       medium
      > 80          good
    """
    if iface_plddt is None:
        return 0.0
    x = iface_plddt

    if x <= 65.0:
        return 0.0
    elif x <= 70.0:
        return _linear_scale(x, 65.0, 70.0, 0.0, 4.0)
    elif x <= 80.0:
        return _linear_scale(x, 70.0, 80.0, 4.0, 8.0)
    else:
        return _linear_scale(x, 80.0, 100.0, 8.0, 10.0)


def _score_cdr_sasa(cdr_sasa: Optional[float]) -> float:
    """
    CDR SASA, higher is better.

    Bands:
      <= 250        fail region
      250 - 300     poor
      300 - 600     medium
      >= 600        good
    """
    if cdr_sasa is None:
        return 0.0
    x = cdr_sasa

    if x <= 250.0:
        return 0.0
    elif x <= 300.0:
        return _linear_scale(x, 250.0, 300.0, 0.0, 4.0)
    elif x <= 600.0:
        return _linear_scale(x, 300.0, 600.0, 4.0, 8.0)
    else:
        return _linear_scale(x, 600.0, 1000.0, 8.0, 10.0)


def _score_netsolp(netsolp: Optional[float]) -> float:
    """
    NetSolP probability in [0,1], higher is better.

    Bands:
      < 0.30        very poor
      0.30 - 0.50   poor
      0.50 - 0.70   medium
      >= 0.70       good
    """
    if netsolp is None:
        return 0.0
    x = _clamp(netsolp, 0.0, 1.0)

    if x <= 0.30:
        return _linear_scale(x, 0.0, 0.30, 0.0, 2.0)
    elif x <= 0.50:
        return _linear_scale(x, 0.30, 0.50, 2.0, 4.0)
    elif x <= 0.70:
        return _linear_scale(x, 0.50, 0.70, 4.0, 8.0)
    else:
        return _linear_scale(x, 0.70, 1.00, 8.0, 10.0)


def _score_novelty(cdr3_identity: Optional[float]) -> float:
    """
    CDR3 identity in percent, lower is better.

    Bands:
      >= 95         worst (copy)
      90 - 95       poor
      70 - 90       medium
      < 70          good
    """
    if cdr3_identity is None:
        return 0.0
    x = cdr3_identity

    if x >= 95.0:
        return 0.0
    elif x >= 90.0:
        return _linear_scale(x, 95.0, 90.0, 0.0, 4.0)
    elif x >= 70.0:
        return _linear_scale(x, 90.0, 70.0, 4.0, 8.0)
    else:
        return _linear_scale(x, 70.0, 0.0, 8.0, 10.0)


def _average(values: List[float]) -> float:
    vals = [v for v in values if v is not None]
    vals = [v for v in vals if v >= 0.0]
    if not vals:
        return 0.0
    return sum(vals) / float(len(vals))


def score_design(raw: RawMetrics) -> FinalScores:
    pm = PerMetricScores()

    pm.delta_g_score = _score_delta_g(raw.delta_g)
    pm.contacts_score = _score_contacts(raw.contacts)
    pm.ipsae_score = _score_ipsae(raw.ipsae)
    pm.dockq_score = _score_dockq(raw.dockq)
    pm.iface_plddt_score = _score_iface_plddt(raw.iface_plddt)
    pm.cdr_sasa_score = _score_cdr_sasa(raw.cdr_sasa)
    pm.netsolp_score = _score_netsolp(raw.netsolp)
    pm.novelty_score = _score_novelty(raw.cdr3_identity)

    # DockQ or ipSAE fallback
    dockq_or_ipsae = pm.dockq_score if raw.dockq is not None else pm.ipsae_score

    binding_values = [
        pm.delta_g_score,
        pm.contacts_score,
        pm.ipsae_score,
        dockq_or_ipsae,
        pm.iface_plddt_score,
        pm.cdr_sasa_score,
    ]
    binding_struct_score = _average(binding_values)

    developability_score = pm.netsolp_score
    novelty_category_score = pm.novelty_score

    final_score_10 = (
        0.60 * binding_struct_score
        + 0.20 * developability_score
        + 0.20 * novelty_category_score
    )
    final_score_100 = final_score_10 * 10.0

    # viability checks
    is_viable = True
    fail_reason = None

    if raw.ipsae is not None and raw.ipsae < 0.60:
        is_viable = False
        fail_reason = "ipSAE < 0.60"

    if is_viable and raw.netsolp is not None and raw.netsolp < 0.50:
        is_viable = False
        fail_reason = "NetSolP < 0.50"

    if is_viable and raw.delta_g is not None and raw.delta_g > -6.0:
        is_viable = False
        fail_reason = "DeltaG > -6 kcal/mol"

    if is_viable and raw.cdr_sasa is not None and raw.cdr_sasa <= 250.0:
        is_viable = False
        fail_reason = "CDR SASA <= 250 Å^2"

    if is_viable and raw.cdr3_identity is not None and raw.cdr3_identity >= 95.0:
        is_viable = False
        fail_reason = "CDR3 identity >= 95 percent"

    category = CategoryScores(
        binding_struct_score=binding_struct_score,
        developability_score=developability_score,
        novelty_category_score=novelty_category_score,
    )

    return FinalScores(
        per_metric=pm,
        category=category,
        final_score_10=final_score_10,
        final_score_100=final_score_100,
        is_viable=is_viable,
        fail_reason=fail_reason,
    )
