"""Generate finite semigroups of partial maps and check Fourier restrictions."""

import argparse
from collections import deque
from itertools import combinations


def closure(generators):
    """Maps use -1 for undefined and index their underlying set from zero."""
    generators = tuple(set(generators))
    seen = set(generators)
    queue = deque(seen)
    while queue:
        current = queue.popleft()
        for step in generators:
            product = tuple(-1 if i < 0 else step[i] for i in current)
            if product not in seen:
                seen.add(product)
                queue.append(product)
    return seen


def fourier_restriction(n, vertices):
    index = {vertex: i for i, vertex in enumerate(vertices)}
    generators = [tuple(index.get((a + g) % n, -1) for a in vertices) for g in range(n)]
    return closure(generators)


def intersections(n, vertices):
    """Independently construct maps via intersections of translates."""
    full = frozenset(vertices)
    domains = {full}
    for shift in range(n):
        allowed = {(a - shift) % n for a in vertices}
        domains.update([domain.intersection(allowed) for domain in domains])
    index = {vertex: i for i, vertex in enumerate(vertices)}
    result = set()
    for domain in domains:
        for shift in range(n):
            if all((a + shift) % n in index for a in domain):
                result.add(tuple(index[(a + shift) % n] if a in domain else -1
                                 for a in vertices))
    return result


def band(s, k):
    vertices = [(a, h) for a in range(s) for h in range(k)]
    index = {vertex: i for i, vertex in enumerate(vertices)}
    generators = [tuple(index.get((a + d, (h + g) % k), -1) for a, h in vertices)
                  for d in range(1 - s, s) for g in range(k)]
    maps = closure(generators)
    fixed_sets = {tuple(i for i, j in enumerate(f) if i == j) for f in maps}
    assert len(maps) == 1 + k * sum(j * j for j in range(1, s + 1))
    assert len(fixed_sets) == 1 + s * (s + 1) // 2
    return len(maps), len(fixed_sets)


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--max-n", type=int, default=8, help="largest cyclic Fourier order")
    args = cli.parse_args()
    for n in range(2, args.max_n + 1):
        for size in range(1, n + 1):
            for vertices in combinations(range(n), size):
                assert fourier_restriction(n, vertices) == intersections(n, vertices)
        punctured = len(fourier_restriction(n, tuple(range(1, n))))
        assert punctured == 2 ** (n - 1) + (n - 1) * (2 ** (n - 2) - 1)
        print("cyclic order", n, "subsets verified", 2 ** n - 1, "punctured size", punctured)
    for s in (2, 3, 4):
        for k in (1, 2, 3):
            print("layers", s, "group order", k, "semigroup and diagonal sizes", band(s, k))
