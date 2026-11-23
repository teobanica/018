"""Exact row-quotient colors for partial Butson matrices, using exponent arrays."""

import argparse
import json
from collections import Counter
from functools import lru_cache
from itertools import combinations, product
from math import lcm


def divide(a, b):
    """Polynomial division by a monic polynomial, least-degree coefficient first."""
    a = list(a)
    quotient = [0] * max(0, len(a) - len(b) + 1)
    for top in range(len(a) - 1, len(b) - 2, -1):
        shift = top - len(b) + 1
        coeff = a[top]
        quotient[shift] = coeff
        for j, value in enumerate(b):
            a[j + shift] -= coeff * value
    return quotient, a[:len(b) - 1]


@lru_cache(None)
def cyclotomic(q):
    coefficients = [-1] + [0] * (q - 1) + [1]
    for d in range(1, q):
        if q % d == 0:
            coefficients, remainder = divide(coefficients, cyclotomic(d))
            assert not any(remainder)
    return tuple(coefficients)


def quotient(a, b, q):
    return tuple((x - y) % q for x, y in zip(a, b))


def vanishes(v, q):
    counts = Counter(v)
    polynomial = [counts[j] for j in range(q)]
    _, remainder = divide(polynomial, cyclotomic(q))
    return not any(remainder)


def analyze(rows, q):
    if q < 1 or not rows or not rows[0]:
        raise ValueError("positive level and a nonempty rectangular array required")
    if any(len(row) != len(rows[0]) or any(type(x) is not int for x in row)
           for row in rows):
        raise ValueError("expected rectangular integer exponent rows")
    normalized = [tuple((x - row[0]) % q for x in row) for row in rows]
    if any(not vanishes(quotient(a, b, q), q) for a, b in combinations(normalized, 2)):
        raise ValueError("input rows are not orthogonal")
    pairs = {}
    for i, j in product(range(len(rows)), repeat=2):
        v = quotient(normalized[i], normalized[j], q)
        pairs.setdefault(v, []).append((i, j))
    vectors = list(pairs)
    neighbors = [set() for _ in vectors]
    for i, j in combinations(range(len(vectors)), 2):
        if not vanishes(quotient(vectors[i], vectors[j], q), q):
            neighbors[i].add(j)
            neighbors[j].add(i)
    commuting = not any(neighbors)
    components = []
    visited = set()
    for start in range(len(vectors)):
        if start in visited:
            continue
        pending = [start]
        component = []
        visited.add(start)
        while pending:
            i = pending.pop()
            component.extend(pairs[vectors[i]])
            for j in sorted(neighbors[i] - visited):
                pending.append(j)
                visited.add(j)
        components.append(sorted(component))
    return {"rows": len(rows), "width": len(rows[0]), "lines": len(pairs),
            "commuting": commuting, "components": components}


def band(s, k, phase_order=7):
    # K is cyclic. The fiber indexed by character c receives phase exp(2*pi*i*c/phase_order).
    t = 2 * s - 1
    q = lcm(t, k, phase_order)
    rows = [[(h * c * (q // k) + a * (c * (q // phase_order) + j * (q // t))) % q
             for c in range(k) for j in range(t)] for a in range(s) for h in range(k)]
    return rows, q


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("input", nargs="?", help='JSON file containing "q" and integer exponent "rows"')
    args = cli.parse_args()
    if args.input:
        with open(args.input, encoding="utf-8") as handle:
            data = json.load(handle)
        print(json.dumps(analyze(data["rows"], data["q"]), indent=2))
    else:
        for s in (2, 3, 4):
            for k in (1, 2, 3):
                rows, q = band(s, k)
                report = analyze(rows, q)
                assert report["commuting"]
                assert len(report["components"]) == report["lines"] == (2 * s - 1) * k
                print("layers", s, "group order", k, "orthogonal quotient lines", report["lines"])
        # Noncommuting 2 x 4 example with row (1,-1,z,-z), z a fifth root.
        report = analyze([[0, 0, 0, 0], [0, 5, 2, 7]], 10)
        assert not report["commuting"] and len(report["components"]) == 2
        print("Noncommuting example:", report["components"])
