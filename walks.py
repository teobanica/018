"""Rational matrix representations separating prescribed quantum coordinate words."""

import argparse
from fractions import Fraction
from itertools import product
import json

from association import association, read_graphs


def gram(n):
    return [[Fraction(i == j) + Fraction(abs(i - j) == 1, 4) for j in range(n)]
            for i in range(n)]


def matrix_product(a, b):
    columns = list(zip(*b))
    return [[sum(x * y for x, y in zip(row, col)) for col in columns] for row in a]


def representation(target_word, neighbors):
    if not target_word or any(i < 0 or i >= len(neighbors) for i in target_word):
        raise ValueError("use a nonempty word in allowed assignments")
    if any(b not in neighbors[a] for a, b in zip(target_word, target_word[1:])):
        raise ValueError("the specified word must be a reduced walk")
    size = len(target_word)
    metric = gram(size)
    # In this basis P_v = E_v B, where E_v selects occurrences of v.
    projections = [[[metric[i][j] if target_word[i] == v else 0 for j in range(size)]
                    for i in range(size)] for v in range(len(neighbors))]
    return metric, projections


def coefficient(target_word, projections, tested_word):
    if not tested_word:
        raise ValueError("first remove the constant by evaluation at the empty map")
    # The last projection sends the dual of the last basis vector to e_last or zero.
    vector = [0] * len(target_word)
    vector[-1] = int(tested_word[-1] == target_word[-1])
    for v in reversed(tested_word[:-1]):
        vector = [sum(x * y for x, y in zip(row, vector)) for row in projections[v]]
    return vector[0]


def verify_representation(metric, projections, neighbors):
    zero = [[0] * len(metric) for _ in metric]
    for v, p in enumerate(projections):
        assert matrix_product(p, p) == p
        assert matrix_product(metric, p) == matrix_product(list(zip(*p)), metric)
        for w in range(v):
            if w not in neighbors[v]:
                assert matrix_product(p, projections[w]) == zero


def run_checks():
    two = [[0, 0], [0, 0]]
    _, neighbors = association(two, two)
    for length in range(1, 5):
        for target in product(range(4), repeat=length):
            if any(b not in neighbors[a] for a, b in zip(target, target[1:])):
                continue
            metric, projections = representation(target, neighbors)
            verify_representation(metric, projections, neighbors)
            for m in range(1, length + 1):
                for tested in product(range(4), repeat=m):
                    expected = Fraction(1, 4 ** (length - 1)) if tested == target else 0
                    assert coefficient(target, projections, tested) == expected
    # Two nonconsecutive incompatible assignments can both occur in a nonzero word.
    three = [[0] * 3 for _ in range(3)]
    vertices, neighbors = association(three, three)
    word = [vertices.index(pair) for pair in [(0, 0), (1, 1), (2, 0)]]
    metric, projections = representation(word, neighbors)
    verify_representation(metric, projections, neighbors)
    assert word[-1] not in neighbors[word[0]]
    assert coefficient(word, projections, word) == Fraction(1, 16)


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("input", nargs="?", help="use the same graph JSON format as association.py")
    cli.add_argument("--walk", help="JSON list of [target,source] assignments, starting at zero")
    args = cli.parse_args()
    if not args.input and not args.walk:
        run_checks()
        print("All projection and word-separation checks passed.")
    else:
        if not args.input or not args.walk:
            cli.error("supply both a graph file and --walk")
        source, target = read_graphs(args.input)
        assignments, neighbors = association(source, target)
        try:
            walk = [assignments.index(tuple(pair)) for pair in json.loads(args.walk)]
            metric, projections = representation(walk, neighbors)
        except (TypeError, ValueError) as err:
            cli.error(str(err))
        verify_representation(metric, projections, neighbors)
        print("Separating coefficient:", coefficient(walk, projections, walk))
        for v in sorted(set(walk)):
            print("Projection for", assignments[v], "in Gram coordinates")
            for row in projections[v]:
                print(*(str(x) for x in row))
