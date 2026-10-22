"""Independent statistical decisions and diagnostic states."""
import math


def decision_fields(q_marginal, q_stratified, analysis_status, balance_status, alpha=0.05):
    if not 0 < alpha < 1:
        raise ValueError('alpha must be in (0, 1)')
    if analysis_status not in {'ok', 'no_randomization', 'constant_statistic', 'not_run'}:
        raise ValueError('unknown analysis_status')
    if balance_status not in {'within_threshold', 'above_threshold', 'undefined'}:
        raise ValueError('unknown balance_status')
    return {
        'marginal_significant': bool(math.isfinite(q_marginal) and q_marginal < alpha),
        'stratified_significant': bool(analysis_status == 'ok' and math.isfinite(q_stratified) and q_stratified < alpha),
        'analysis_status': analysis_status,
        'balance_status': balance_status,
    }
