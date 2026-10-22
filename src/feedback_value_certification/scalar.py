"""Exact scalar value evaluation on rational cuts within reachable envelopes."""

from fractions import Fraction as F
from bisect import bisect_left


def _require(condition, message):
    if not condition:
        raise ValueError(message)


class EvaluationLimitError(RuntimeError):
    """A requested partition exceeds the declared construction limit."""


def q(v):
    if isinstance(v, float):
        raise TypeError("Use exact rational strings, not floats")
    return F(v)


class Model:

    def __init__(self, raw, policy):
        self.name, self.policy = (raw["name"], policy)
        self.lo, self.hi = map(q, raw["domain"])
        self.lam = q(raw["contraction"])
        self.x0 = q(raw["x0"])
        self.v = list(map(q, raw["increments"]))
        self.S = len(self.v)
        self.init = list(map(q, raw["initial_state_law"]))
        self.thresholds = list(map(q, raw["thresholds"]))
        self.point_tags = list(raw["point_tags"])
        self.open_tags = list(raw["open_tags"])
        self.mu = [list(map(q, row)) for row in raw["policies"][policy]]
        self.P = [[[q(v) for v in row] for row in actions] for actions in raw["response"]]
        self.reward = [[[q(v) for v in row] for row in actions] for actions in raw["reward"]]
        _require(
            self.lo < self.hi and 0 < self.lam < 1 and (self.lo <= self.x0 <= self.hi),
            "Invalid scalar model or inconsistent value partition",
        )
        _require(
            self.thresholds == sorted(set(self.thresholds)),
            "Invalid scalar model or inconsistent value partition",
        )
        _require(
            all((self.lo < x < self.hi for x in self.thresholds)),
            "Invalid scalar model or inconsistent value partition",
        )
        _require(
            len(self.point_tags) == len(self.thresholds),
            "Invalid scalar model or inconsistent value partition",
        )
        _require(
            len(self.open_tags) == len(self.thresholds) + 1,
            "Invalid scalar model or inconsistent value partition",
        )
        _require(
            len(self.init) == self.S and sum(self.init) == 1 and (min(self.init) >= 0),
            "Invalid scalar model or inconsistent value partition",
        )
        _require(
            all((sum(row) == 1 and min(row) >= 0 for row in self.mu)),
            "Invalid scalar model or inconsistent value partition",
        )
        A = len(self.mu[0])
        _require(
            all((len(row) == A for row in self.mu)),
            "Invalid scalar model or inconsistent value partition",
        )
        _require(
            len(self.P) == len(self.reward) == self.S,
            "Invalid scalar model or inconsistent value partition",
        )
        for z in range(self.S):
            _require(
                len(self.P[z]) == len(self.reward[z]) == A,
                "Invalid scalar model or inconsistent value partition",
            )
            for a in range(A):
                _require(
                    len(self.P[z][a]) == len(self.reward[z][a]) == self.S,
                    "Invalid scalar model or inconsistent value partition",
                )
                _require(
                    sum(self.P[z][a]) == 1 and min(self.P[z][a]) >= 0,
                    "Invalid scalar model or inconsistent value partition",
                )
                _require(
                    max(map(abs, self.reward[z][a])) <= 1,
                    "Invalid scalar model or inconsistent value partition",
                )
            _require(
                self.lo <= self.image(self.lo, z) <= self.image(self.hi, z) <= self.hi,
                "Invalid scalar model or inconsistent value partition",
            )
        tags = set(self.point_tags + self.open_tags)
        _require(
            tags <= set(range(len(self.mu))), "Invalid scalar model or inconsistent value partition"
        )
        self.atoms = []
        cuts = [self.lo] + self.thresholds + [self.hi]
        for i, x in enumerate(cuts):
            self.atoms.append((x, x, self.tag(x), False))
            if i + 1 < len(cuts):
                self.atoms.append((x, cuts[i + 1], self.open_tags[i], True))
        self.K, self.J, self.R = ([], [], [])
        for dist in self.mu:
            kk = []
            jj = []
            rr = []
            for z in range(self.S):
                kr = []
                jr = []
                for zp in range(self.S):
                    kr.append(sum((dist[a] * self.P[z][a][zp] for a in range(A))))
                    jr.append(
                        sum((dist[a] * self.P[z][a][zp] * self.reward[z][a][zp] for a in range(A)))
                    )
                _require(sum(kr) == 1, "Invalid scalar model or inconsistent value partition")
                kk.append(kr)
                jj.append(jr)
                rr.append(sum(jr))
            self.K.append(kk)
            self.J.append(jj)
            self.R.append(rr)

    def image(self, x, z):
        return self.lam * x + self.v[z]

    def tag(self, x):
        i = bisect_left(self.thresholds, x)
        return (
            self.point_tags[i]
            if i < len(self.thresholds) and self.thresholds[i] == x
            else self.open_tags[i]
        )

    def intersections(self, interval):
        """Actual point/open atoms intersected with a CLOSED envelope.

        Return closed hulls, but only after checking actual intersection is
        nonempty. A zero-width envelope inside an open atom remains nonempty.
        """
        l, r = interval
        for u, v, tag, strict in self.atoms:
            a, b = (max(l, u), min(r, v))
            if a > b:
                continue
            if strict and a == b and (not u < a < v):
                continue
            yield (a, b, tag)


class Function:

    def __init__(self, cuts, points, opens):
        _require(
            cuts and len(points) == len(cuts) and (len(opens) == len(cuts) - 1),
            "Invalid scalar model or inconsistent value partition",
        )
        self.cuts, self.points, self.opens = (cuts, points, opens)

    def at(self, x):
        i = bisect_left(self.cuts, x)
        if i < len(self.cuts) and self.cuts[i] == x:
            return self.points[i]
        _require(0 < i < len(self.cuts), ("continuation outside envelope", str(x)))
        return self.opens[i - 1]

    def coalesce(self):
        if len(self.cuts) <= 2:
            return self
        keep = (
            [0]
            + [
                i
                for i in range(1, len(self.cuts) - 1)
                if not self.points[i] == self.opens[i - 1] == self.opens[i]
            ]
            + [len(self.cuts) - 1]
        )
        return Function(
            [self.cuts[i] for i in keep],
            [self.points[i] for i in keep],
            [self.opens[i] for i in keep[:-1]],
        )


def envelopes(model, T):
    levels = [[(model.x0, model.x0) if mass else None for mass in model.init]]
    for t in range(T):
        nxt = [None] * model.S
        for z, interval in enumerate(levels[-1]):
            if interval is None:
                continue
            for lo, hi, tag in model.intersections(interval):
                a, b = (model.image(lo, z), model.image(hi, z))
                for zp, prob in enumerate(model.K[tag][z]):
                    if not prob:
                        continue
                    old = nxt[zp]
                    nxt[zp] = (a, b) if old is None else (min(old[0], a), max(old[1], b))
        _require(
            any((v is not None for v in nxt)),
            "Invalid scalar model or inconsistent value partition",
        )
        levels.append(nxt)
    return levels


def evaluate_exact(model, horizon, *, max_cuts=12000, max_atoms=200000):
    """Compute an exact rational finite-horizon value with actual threshold tags.

    Return a JSON-compatible value, singleton interval, and construction counts.
    """
    if isinstance(horizon, bool) or not isinstance(horizon, int) or horizon < 0:
        raise ValueError("Horizon must be a nonnegative integer")
    if any((isinstance(x, bool) or not isinstance(x, int) or x < 1 for x in (max_cuts, max_atoms))):
        raise ValueError("Construction limits must be positive integers")
    T = horizon
    guards = {
        "max_cuts_per_state_layer": max_cuts,
        "max_evaluated_atoms_per_method_call": max_atoms,
    }
    stats = {
        "constructed_atoms": 0,
        "evaluated_point_atoms": 0,
        "evaluated_open_atoms": 0,
        "retained_atoms": 0,
        "removed_breakpoints": 0,
        "max_cuts_state_layer": 0,
        "positive_kernel_terms": 0,
        "envelope_nonempty_state_layers": 0,
        "method": "symbolic_reachable",
    }
    levels = envelopes(model, T)
    stats["envelope_nonempty_state_layers"] = sum((x is not None for row in levels for x in row))
    nxt = [
        (
            None
            if v is None
            else Function(
                [v[0]] if v[0] == v[1] else list(v),
                [F(0)] if v[0] == v[1] else [F(0), F(0)],
                [] if v[0] == v[1] else [F(0)],
            )
        )
        for v in levels[T]
    ]
    for t in range(T - 1, -1, -1):
        current = []
        for z, interval in enumerate(levels[t]):
            if interval is None:
                current.append(None)
                continue
            lo, hi = interval
            possible_tags = {tag for a, b, tag in model.intersections(interval)}
            cuts = {lo, hi}
            cuts.update((x for x in model.thresholds if lo < x < hi))
            for zp in range(model.S):
                if any((model.K[tag][z][zp] for tag in possible_tags)):
                    _require(
                        nxt[zp] is not None, "Invalid scalar model or inconsistent value partition"
                    )
                    for y in nxt[zp].cuts:
                        x = (y - model.v[z]) / model.lam
                        if lo < x < hi:
                            cuts.add(x)
            cuts = sorted(cuts)
            if len(cuts) > guards["max_cuts_per_state_layer"]:
                raise EvaluationLimitError("CUT_CAP at time%d,state%d,cuts%d" % (t, z, len(cuts)))
            stats["max_cuts_state_layer"] = max(stats["max_cuts_state_layer"], len(cuts))
            points = []
            opens = []

            def evaluate(x):
                tag = model.tag(x)
                y = model.image(x, z)
                value = model.R[tag][z]
                for zp, prob in enumerate(model.K[tag][z]):
                    if prob:
                        _require(
                            nxt[zp] is not None,
                            "Invalid scalar model or inconsistent value partition",
                        )
                        value += prob * nxt[zp].at(y)
                        stats["positive_kernel_terms"] += 1
                return value

            stats["constructed_atoms"] += 2 * len(cuts) - 1
            if stats["constructed_atoms"] > guards["max_evaluated_atoms_per_method_call"]:
                raise EvaluationLimitError("TOTAL_ATOM_CAP at time%d" % t)
            for i, x in enumerate(cuts):
                points.append(evaluate(x))
                stats["evaluated_point_atoms"] += 1
                if i + 1 < len(cuts):
                    opens.append(evaluate((x + cuts[i + 1]) / 2))
                    stats["evaluated_open_atoms"] += 1
            before = len(cuts)
            fn = Function(cuts, points, opens).coalesce()
            stats["removed_breakpoints"] += before - len(fn.cuts)
            stats["retained_atoms"] += 2 * len(fn.cuts) - 1
            current.append(fn)
        nxt = current
    value = sum((mass * nxt[z].at(model.x0) for z, mass in enumerate(model.init) if mass))
    return {
        "status": "COMPLETE_EXACT",
        "value": str(value),
        "interval": [str(value), str(value)],
        "width": "0",
        "stats": stats,
    }
