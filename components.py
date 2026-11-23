"""Check adjacency commutation and decomposition over connected components."""

from itertools import combinations, permutations

from diagonals import difference, multiply, projector


def components(neighbors):
    remaining = set(range(len(neighbors)))
    result = []
    while remaining:
        seen = set()
        todo = [min(remaining)]
        while todo:
            vertex = todo.pop()
            if vertex not in seen:
                seen.add(vertex)
                todo.extend(neighbors[vertex] - seen)
        remaining -= seen
        result.append(frozenset(seen))
    return result


def partial_maps(n):
    for rank in range(n + 1):
        for domain in combinations(range(n), rank):
            for image in permutations(range(n), rank):
                yield dict(zip(domain, image))


def commutes(neighbors, mapping):
    # Compare the supports of the two columns at a in Du and uD.
    for a, adjacent in enumerate(neighbors):
        left = neighbors[mapping[a]] if a in mapping else set()
        right = {mapping[b] for b in adjacent if b in mapping}
        if left != right:
            return False
    return True


def is_component_map(neighbors, parts, mapping):
    for source in parts:
        domain = source.intersection(mapping)
        if not domain:
            continue
        if domain != source or {mapping[a] for a in source} not in parts:
            return False
        for a, b in combinations(source, 2):
            if (a in neighbors[b]) != (mapping[a] in neighbors[mapping[b]]):
                return False
    return True


def check_scalar_graphs(n):
    edges = list(combinations(range(n), 2))
    maps = list(partial_maps(n))
    accepted = 0
    for mask in range(1 << len(edges)):
        neighbors = [set() for _ in range(n)]
        for bit, (a, b) in enumerate(edges):
            if mask & (1 << bit):
                neighbors[a].add(b)
                neighbors[b].add(a)
        parts = components(neighbors)
        for mapping in maps:
            works = commutes(neighbors, mapping)
            assert works == is_component_map(neighbors, parts, mapping)
            accepted += works
    return len(maps) * (1 << len(edges)), accepted


ZERO = [[0, 0], [0, 0]]
UNIT = [[1, 0], [0, 1]]


def matrix_sum(matrices):
    result = [[0, 0], [0, 0]]
    for matrix in matrices:
        for i in range(2):
            for j in range(2):
                result[i][j] += matrix[i][j]
    return result


def verify_array(neighbors, entries):
    """Check the defining relations and recover every block support two ways."""
    n = len(neighbors)
    def at(i, a):
        return entries.get((i, a), ZERO)

    for value in entries.values():
        assert value == list(map(list, zip(*value)))
        assert multiply(value, value) == value
    for i in range(n):
        row = [value for (j, a), value in entries.items() if j == i]
        col = [value for (a, j), value in entries.items() if j == i]
        for line in (row, col):
            assert all(multiply(left, right) == ZERO for left, right in combinations(line, 2))
        for a in range(n):
            assert matrix_sum(at(j, a) for j in neighbors[i]) == matrix_sum(
                at(i, b) for b in neighbors[a])

    supports = []
    parts = components(neighbors)
    for target in parts:
        row = []
        for source in parts:
            sums = [matrix_sum(at(i, a) for a in source) for i in target]
            sums.extend(matrix_sum(at(i, a) for i in target) for a in source)
            assert all(value == sums[0] for value in sums)
            row.append(sums[0])
        supports.append(row)
    for line in [*supports, *zip(*supports)]:
        assert all(multiply(left, right) == ZERO for left, right in combinations(line, 2))
    return supports


def check_quantum_blocks():
    # Two disjoint edges exchange on the complementary line. A third,
    # isolated component acts on a line that does not commute with either.
    neighbors = [{1}, {0}, {3}, {2}, set()]
    p = projector([1, 0])
    q = projector([1, 1])
    other = difference(UNIT, p)
    entries = {(4, 4): q}
    for a in range(4):
        entries[a, a] = p
        entries[(a + 2) % 4, a] = other
    supports = verify_array(neighbors, entries)
    assert supports == [[p, other, ZERO], [other, p, ZERO], [ZERO, ZERO, q]]
    assert multiply(p, q) != multiply(q, p)


def check_decorated_cycle(n):
    # Periodic model of the decorated path. At displacement h the same
    # four-point magic matrix acts at each site. Test each translation fiber.
    neighbors = [set() for _ in range(5 * n)]
    for j in range(n):
        root = 5 * j
        neighbors[root].update({5 * ((j - 1) % n), 5 * ((j + 1) % n)})
        for leaf in range(1, 5):
            neighbors[root].add(root + leaf)
            neighbors[root + leaf].add(root)
    p = projector([1, 0])
    q = projector([1, 1])
    for h in range(n):
        entries = {}
        for j in range(n):
            source, target = 5 * j, 5 * ((j + h) % n)
            entries[target, source] = UNIT
            for start, line in [(1, p), (3, q)]:
                for offset in (0, 1):
                    a = source + start + offset
                    entries[target + start + offset, a] = line
                    entries[target + start + 1 - offset, a] = difference(UNIT, line)
        assert verify_array(neighbors, entries) == [[UNIT]]


if __name__ == "__main__":
    for order in range(1, 5):
        tested, passing = check_scalar_graphs(order)
        print(f"{order} vertices: {tested} graph/map pairs checked, {passing} commute")
    check_quantum_blocks()
    check_decorated_cycle(5)
    print("Recovered all matrix-valued block supports and verified the five-site periodic model.")
