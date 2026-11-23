"""Cliques and walks for partial isomorphisms of finite colored directed graphs."""

import argparse
from itertools import combinations, permutations, product
import json


def check_square(matrix):
    if any(len(row) != len(matrix) for row in matrix):
        raise ValueError("each graph must be a square color matrix")


def read_graphs(filename):
    with open(filename, encoding="utf-8") as handle:
        data = json.load(handle)
    source = data["source"]
    return source, data.get("target", source)


def association(source, target):
    """Return assignments (target, source) and adjacency sets in this order."""
    check_square(source)
    check_square(target)
    assignments = [(i, a) for i in range(len(target)) for a in range(len(source))
                   if target[i][i] == source[a][a]]
    neighbors = [set() for v in assignments]
    for left, right in combinations(range(len(assignments)), 2):
        i, a = assignments[left]
        j, b = assignments[right]
        if i == j or a == b:
            continue
        if (target[i][j], target[j][i]) == (source[a][b], source[b][a]):
            neighbors[left].add(right)
            neighbors[right].add(left)
    return assignments, neighbors


def clique_counts(neighbors, max_rank):
    """A clique is generated once, by adjoining vertices in increasing order."""
    totals = [1] + [0] * max_rank

    def extend(candidates, size):
        for pos, vertex in enumerate(candidates):
            totals[size + 1] += 1
            suffix = [other for other in candidates[pos + 1:] if other in neighbors[vertex]]
            if suffix:
                extend(suffix, size + 1)

    extend(list(range(len(neighbors))), 0)
    return totals


def walk_counts(neighbors, max_length):
    counts = [1]  # The unit has degree zero, including for an empty association graph.
    ending = [1] * len(neighbors)
    for length in range(1, max_length + 1):
        counts.append(sum(ending))
        ending = [sum(ending[j] for j in adjacent) for adjacent in neighbors]
    return counts


def direct_partial_maps(source, target):
    """Independent verification from restrictions of labeled partial bijections."""
    totals = []
    for rank in range(min(len(source), len(target)) + 1):
        count = 0
        for domain in combinations(range(len(source)), rank):
            for image in permutations(range(len(target)), rank):
                if all(source[a][b] == target[i][j]
                       for a, i in zip(domain, image) for b, j in zip(domain, image)):
                    count += 1
        totals.append(count)
    return totals


def example_checks():
    small = [[[bits[0], bits[1]], [bits[2], bits[3]]] for bits in product((0, 1), repeat=4)]
    for source, target in product(small, repeat=2):
        assignments, neighbors = association(source, target)
        assert clique_counts(neighbors, 2) == direct_partial_maps(source, target)
        if source == target:
            if source[0][0] != source[1][1]:
                expected = 4
            elif source[0][1] != source[1][0]:
                expected = 6
            else:
                expected = 7
            assert sum(clique_counts(neighbors, 2)) == expected
    edge = [[0, 0, 0], [0, 0, 1], [0, 1, 0]]
    _, neighbors = association(edge, edge)
    assert walk_counts(neighbors, 2) == [1, 9, 20]
    assert clique_counts(neighbors, 3) == [1, 9, 10, 2]
    for n in range(3, 8):
        for directed in (True, False):
            matrix = [[int((j-i) % n == 1 or (not directed and (i-j) % n == 1))
                       for j in range(n)] for i in range(n)]
            _, neighbors = association(matrix, matrix)
            degree = (2 if directed else 4) + (n - 3) ** 2
            assert walk_counts(neighbors, 5) == [1] + [n*n*degree**j for j in range(5)]
    # Different numbers of vertices, colored loops, and no possible assignment.
    for source, target in [(small[5], edge), (edge, small[2]), ([[2]], [[3]]), ([], edge)]:
        _, neighbors = association(source, target)
        assert clique_counts(neighbors, min(len(source), len(target))) == direct_partial_maps(source, target)


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("input", nargs="?", help="JSON with a source color matrix and optional target")
    cli.add_argument("--degree", type=int, default=5, help="maximum word length")
    args = cli.parse_args()
    if args.degree < 0:
        cli.error("degree must be nonnegative")
    if args.input:
        x, y = read_graphs(args.input)
    else:
        example_checks()
        print("All checks passed. Reporting the single-edge graph on three points.")
        x = y = [[0, 0, 0], [0, 0, 1], [0, 1, 0]]
    vertices, graph = association(x, y)
    print("Allowed assignments (target, source):", vertices)
    print("Classical counts by rank:", clique_counts(graph, min(len(x), len(y))))
    print("Quantum basis counts by length:", walk_counts(graph, args.degree))
