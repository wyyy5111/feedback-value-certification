"""Exact rational certificates for projected group-word probabilities.

Equal projected digits are grouped before the finite-state survival check.
"""

from fractions import Fraction as F
from math import lcm


def exact(x):
    if isinstance(x, float):
        raise TypeError("Use exact rational inputs")
    return F(x)


def projected_certificate(kernels, projections, contraction):
    """Derive complete collision groups from exact boundary digits.

    projections contains one obs->affine-digit map per boundary. Callers are
    responsible for deriving the digits from their actual affine update.
    A residue check supplies one sufficient rational word-separation rule.
    """
    if not kernels:
        raise ValueError("Provide at least one control tag")
    factor = exact(contraction)
    if not 0 < factor < 1:
        raise ValueError("Contraction must be between zero and one")
    observations = tuple(next(iter(kernels.values())))
    if not projections:
        raise ValueError("No boundary projections")
    categories, geometry = [], []
    for projection in projections:
        if set(projection) != set(observations):
            raise ValueError("Projection observation alphabet differs")
        groups = {}
        for obs in observations:
            groups.setdefault(exact(projection[obs]), []).append(obs)
        digits = list(groups)
        reference = digits[0]
        denominator = lcm(*[(x - reference).denominator for x in digits])
        integers = [int((x - reference) * denominator) for x in digits]
        residues = [x % factor.denominator for x in integers]
        for group in groups.values():
            if group not in categories:
                categories.append(group)
        geometry.append(
            {
                "distinct_projected_digits": [str(x) for x in digits],
                "collision_groups": list(groups.values()),
                "centering_denominator_M": denominator,
                "centered_integer_digits": integers,
                "distinct_residues_mod_contraction_denominator": len(set(residues))
                == len(residues),
                "shift_separation_for_N_symbols": "1/(M*q^(N-1))",
            }
        )
    result = certificate(kernels, categories)
    result["projection_geometry"] = geometry
    result["contraction"] = str(factor)
    if not all(g["distinct_residues_mod_contraction_denominator"] for g in geometry):
        result["survival_certificate_status"] = result["status"]
        result["status"] = "NO_CERTIFICATE_RATIONAL_WORD_SEPARATION"
    return result


def certificate(kernels, categories):
    """kernels[tag][observation][source_state][target_state].

    Output categories may overlap. A fixed projected-symbol word uses one
    category at each step. We upper bound it by allowing category AND tag to
    depend on the full current finite state; this extra power is conservative.
    A category must include ALL raw symbols with the same boundary projection.
    """
    tags = list(kernels)
    if not tags:
        raise ValueError("No control tags")
    observations = tuple(kernels[tags[0]])
    if not observations:
        raise ValueError("No observations")
    n = len(kernels[tags[0]][observations[0]])
    if n < 1:
        raise ValueError("No finite states")
    data = {}
    for tag in tags:
        if set(kernels[tag]) != set(observations):
            raise ValueError("Observation alphabets differ")
        data[tag] = {}
        for obs in observations:
            matrix = kernels[tag][obs]
            if len(matrix) != n or any(len(row) != n for row in matrix):
                raise ValueError("Matrix dimensions differ")
            mat = tuple(tuple(exact(x) for x in row) for row in matrix)
            if any(x < 0 for row in mat for x in row):
                raise ValueError("Negative probability")
            data[tag][obs] = mat
        for z in range(n):
            if sum(data[tag][obs][z][v] for obs in observations for v in range(n)) != 1:
                raise ValueError("Controlled kernel row not normalized")
    cats = [tuple(cat) for cat in categories]
    if not cats or any(
        not cat or len(set(cat)) != len(cat) or not set(cat) <= set(observations) for cat in cats
    ):
        raise ValueError("Invalid output categories")
    actions = []
    for tag in tags:
        for cat in cats:
            mat = tuple(
                tuple(sum(data[tag][obs][z][v] for obs in cat) for v in range(n)) for z in range(n)
            )
            actions.append((tag, cat, mat))
    surviving = set(range(n))
    removals = []
    while surviving:
        keep = {
            z
            for z in surviving
            if any(sum(mat[z][v] for v in surviving) == 1 for _, _, mat in actions)
        }
        removed = sorted(surviving - keep)
        if not removed:
            witness = {}
            for z in sorted(surviving):
                tag, cat, _ = next(a for a in actions if sum(a[2][z][v] for v in surviving) == 1)
                witness[str(z)] = {"control_tag": tag, "category": list(cat)}
            return {
                "status": "NO_CERTIFICATE_CLOSED_SURVIVAL_CLASS",
                "surviving_states": sorted(surviving),
                "state_adaptive_survival_witness": witness,
                "removed_layers": removals,
                "not_a_counterexample_to_true_boundary_anti_concentration": True,
            }
        removals.append(removed)
        surviving = keep
    # Removal rank implies absorption within at most m steps has positive
    # probability for EVERY state/tag/category selector. Compute the sharper
    # exact survival bound rather than using a product of minimum probabilities.
    m = len(removals)
    values = (F(1),) * n
    for _ in range(m):
        values = tuple(
            max(sum(mat[z][v] * values[v] for v in range(n)) for _, _, mat in actions)
            for z in range(n)
        )
    q = max(values)
    if not 0 <= q < 1:
        raise ArithmeticError("Closed-set certificate contradicts block survival")
    return {
        "status": "CERTIFIED_UNIFORM_BLOCK_WORD_DECAY",
        "finite_states": n,
        "block_length": m,
        "block_survival_upper": str(q),
        "statewise_block_upper": [str(x) for x in values],
        "removed_layers": removals,
        "word_probability_bound": "q^floor(k/m) for k future emitted symbols; q^floor((N-1)/m) for N-symbol words whose first symbol is fixed",
        "control_and_category_are_allowed_state_adaptive_for_upper_bound": True,
    }
