"""Count partial automorphisms of directed and undirected cycles by rank."""

import argparse
from collections import Counter
from itertools import combinations, permutations
from math import comb, factorial, prod


def partitions(total, minimum=1):
    if total == 0:
        yield ()
    for first in range(minimum, total + 1):
        for tail in partitions(total - first, first):
            yield (first,) + tail


def formula(n, directed):
    counts = [1] + [0] * (n - 1) + [n if directed else 2 * n]
    for size in range(1, n):
        for runs in partitions(size):
            parts = len(runs)
            if parts > n - size:
                continue
            mult = Counter(runs)
            reorderings = prod(factorial(freq) for freq in mult.values())
            marked = n * factorial(parts - 1) * comb(n - size - 1, parts - 1)
            subsets, remainder = divmod(marked, reorderings)
            assert remainder == 0
            reversals = 1 if directed else 2 ** (parts - mult[1])
            counts[size] += subsets ** 2 * reorderings * reversals
    return counts


def component_counts(n):
    """Count subsets by the number of cyclic components, including both constant subsets."""
    return [1] + [2 * comb(n, 2 * p) + int(p == 1) for p in range(1, n // 2 + 1)]


def check_component_counts(n):
    observed = [0] * (1 + n // 2)
    for mask in range(1 << n):
        starts = sum(bool(mask & (1 << i)) and not (mask & (1 << ((i - 1) % n)))
                     for i in range(n))
        count = 1 if mask == (1 << n) - 1 else starts
        observed[count] += 1
    assert observed == component_counts(n)


def brute_force(n, directed):
    def adjacent(i, j):
        return (j - i) % n == 1 or (not directed and (i - j) % n == 1)

    def signature(vertices):
        return tuple(adjacent(i, j) for i in vertices for j in vertices)

    counts = []
    for rank in range(n + 1):
        # The order on the domain gives each bijection a unique encoding.
        domains = Counter(signature(v) for v in combinations(range(n), rank))
        images = Counter(signature(v) for v in permutations(range(n), rank))
        counts.append(sum(freq * images[key] for key, freq in domains.items()))
    return counts


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--max-n", type=int, default=8)
    cli.add_argument("--brute-through", type=int, default=7)
    args = cli.parse_args()
    if args.max_n < 3:
        cli.error("cycles have at least three vertices")
    for order in range(3, args.max_n + 1):
        if order <= args.brute_through:
            check_component_counts(order)
        for oriented in (True, False):
            result = formula(order, oriented)
            if order <= args.brute_through:
                assert result == brute_force(order, oriented)
            print(order, "directed" if oriented else "undirected", result, sum(result))
