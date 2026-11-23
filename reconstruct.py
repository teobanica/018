"""Recover a compatibility graph from equal sums of pairs of characters."""
from collections import Counter
from itertools import combinations, permutations

from association import association, clique_counts, direct_partial_maps


def character_points(graph):
    """Encode clique indicators in base three. Adding two causes no carries."""
    weights = [3**v for v in range(len(graph))]
    points = []

    def extend(value, candidates):
        points.append(value)
        for j, v in enumerate(candidates):
            others = [w for w in candidates[j + 1:] if w in graph[v]]
            extend(value + weights[v], others)

    extend(0, list(range(len(graph))))
    return sorted(points, reverse=True)


def recovered_graphs(points):
    """Use only pair sums. Try every possible base of minimum degree."""
    multiplicity = Counter(a + b for a, b in combinations(points, 2))
    links = {a: set() for a in points}
    for a, b in combinations(points, 2):
        if multiplicity[a + b] == 1:
            links[a].add(b)
            links[b].add(a)
    degree = min(map(len, links.values()))
    for adjacent in links.values():
        if len(adjacent) == degree:
            ordered = sorted(adjacent)
            yield [{j for j, b in enumerate(ordered) if a != b and b not in links[a]}
                   for a in ordered]


def isomorphic(left, right):
    """Check the recovered graph by independent vertex relabeling."""
    if len(left) != len(right):
        return False
    for order in permutations(range(len(left))):
        if all({order[w] for w in neighbors} == right[order[v]]
               for v, neighbors in enumerate(left)):
            return True
    return False


def check_small_graphs():
    checked = 0
    for n in range(6):
        pairs = list(combinations(range(n), 2))
        for bits in range(1 << len(pairs)):
            graph = [set() for v in range(n)]
            for j, (a, b) in enumerate(pairs):
                if bits & (1 << j):
                    graph[a].add(b)
                    graph[b].add(a)
            for result in recovered_graphs(character_points(graph)):
                assert isomorphic(graph, result), (n, bits)
            checked += 1
    print(f"Recovered all {checked} labeled graphs on at most five vertices.")


def check_squares():
    for oriented in (True, False):
        edges = (1,) if oriented else (1, 3)
        matrix = [[int((b - a) % 4 in edges) for b in range(4)] for a in range(4)]
        _, graph = association(matrix, matrix)
        expected = direct_partial_maps(matrix, matrix)
        for recovered in recovered_graphs(character_points(graph)):
            counts = clique_counts(recovered, len(recovered))
            assert counts == expected + [0] * (len(counts) - len(expected))
        label = "Oriented" if oriented else "Unoriented"
        print(f"{label} square, recovered rank counts: {expected}")


if __name__ == "__main__":
    check_small_graphs()
    check_squares()
