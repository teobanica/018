"""Exact fixed-rank fixed-point laws, checked by enumerating partial bijections."""

import argparse
from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations
from math import comb, factorial


def rank_law(n, k, length):
    """Coefficients in z of the generating polynomial in the paper."""
    law = [Fraction(0) for _ in range(min(k, length) + 1)]
    for a in range(len(law)):
        weight = Fraction(comb(k, a) * comb(length, a), factorial(a) * comb(n, a) ** 2)
        for count in range(a + 1):
            law[count] += (-1) ** (a - count) * comb(a, count) * weight
    return law


def direct_laws(n, k):
    counts = [Counter() for _ in range(n + 1)]
    for domain in combinations(range(n), k):
        for image in permutations(range(n), k):
            fixed = {i for i, j in zip(domain, image) if i == j}
            value = 0
            counts[0][0] += 1
            for position in range(n):
                value += position in fixed
                counts[position + 1][value] += 1
    return counts


def partitions(degree):
    """Restricted-growth strings, one per partition of a nonempty set."""
    def extend(blocks):
        if len(blocks) == degree:
            yield tuple(blocks)
            return
        for label in range(max(blocks) + 2):
            yield from extend(blocks + [label])
    yield from extend([0])


def noncrossing(labels):
    for a, b, c, d in combinations(range(len(labels)), 4):
        if labels[a] == labels[c] and labels[b] == labels[d] and labels[a] != labels[b]:
            return False
    return True


def join_size(left, right):
    parent = list(range(len(left)))

    def root(i):
        while parent[i] != i:
            i = parent[i]
        return i

    for pattern in (left, right):
        representatives = {}
        for i, label in enumerate(pattern):
            first = representatives.setdefault(label, i)
            parent[root(i)] = root(first)
    return len({root(i) for i in range(len(left))})


def inverse(matrix):
    size = len(matrix)
    rows = [[Fraction(v) for v in row] + [Fraction(i == j) for j in range(size)]
            for i, row in enumerate(matrix)]
    for j in range(size):
        pivot = next(i for i in range(j, size) if rows[i][j])
        rows[j], rows[pivot] = rows[pivot], rows[j]
        scale = rows[j][j]
        rows[j] = [v / scale for v in rows[j]]
        for i in range(size):
            if i != j and rows[i][j]:
                scale = rows[i][j]
                rows[i] = [x - scale * y for x, y in zip(rows[i], rows[j])]
    return [row[size:] for row in rows]


def finite_moment(n, k, length, degree, quantum):
    from diagonals import multiply

    labels = [p for p in partitions(degree) if not quantum or noncrossing(p)]
    joins = [[join_size(p, q) for q in labels] for p in labels]

    def gram(size):
        return [[size ** blocks for blocks in row] for row in joins]

    w = inverse(gram(n))
    left = multiply(w, gram(k))
    right = multiply(w, gram(length))
    return sum(value * right[j][i] for i, row in enumerate(left) for j, value in enumerate(row))


def check_moments():
    n, k, length = 6, 4, 3
    law = rank_law(n, k, length)
    for d in range(1, 5):
        exact = sum(prob * count ** d for count, prob in enumerate(law))
        assert finite_moment(n, k, length, d, False) == exact
        if d < 4:
            assert finite_moment(n, k, length, d, True) == exact
    assert finite_moment(4, 4, 4, 4, True) == 14
    assert finite_moment(4, 4, 4, 4, False) == 15


def verify(limit):
    tests = 0
    for n in range(limit + 1):
        for k in range(n + 1):
            counts = direct_laws(n, k)
            total = comb(n, k) ** 2 * factorial(k)
            for length, observed in enumerate(counts):
                predicted = rank_law(n, k, length)
                assert sum(observed.values()) == total
                assert sum(predicted) == 1
                assert all(observed[j] == p * total for j, p in enumerate(predicted))
                tests += 1
    return tests


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--moment", nargs=4, type=int, metavar=("N", "K", "L", "D"))
    mode.add_argument("--law", nargs=3, type=int, metavar=("N", "K", "L"))
    args = parser.parse_args()
    if args.moment is not None:
        n, k, length, d = args.moment
        if not (1 <= d <= 5 and n >= k >= max(4, d) and 0 <= length <= n):
            parser.error("require 1 <= D <= 5, N >= K >= max(4,D), and 0 <= L <= N")
        print("Classical moment:", finite_moment(n, k, length, d, False))
        print("Quantum moment:", finite_moment(n, k, length, d, True))
    elif args.law is None:
        check_moments()
        print("Verified the finite Gram trace against classical moments and the fourth Haar moments.")
        print(f"Verified {verify(6)} rank and truncation cases for sets of up to six points.")
    else:
        n, k, length = args.law
        if not (0 <= k <= n and 0 <= length <= n):
            parser.error("require 0 <= K,L <= N")
        for count, probability in enumerate(rank_law(n, k, length)):
            print(count, probability)
