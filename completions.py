"""Minimum magic completions over M_2 and the 18-vertex coloring example.

A positive integer and its negative denote complementary lines; different
absolute values specify different orthogonal bases. Zero and ONE denote the
zero and identity projections. No approximation of angles is needed.
"""

from collections import Counter
from functools import lru_cache
from itertools import combinations, permutations, product

ONE = "I"


def orthogonal(left, right):
    if left == 0 or right == 0:
        return True
    return left != ONE and right != ONE and left == -right


def submagic_line(entries):
    return all(orthogonal(a, b) for a, b in combinations(entries, 2))


def deficit(entries):
    nonzero = [x for x in entries if x != 0]
    if not nonzero:
        return ONE
    if len(nonzero) == 1 and nonzero[0] != ONE:
        return -nonzero[0]
    return 0


def minimum_size(corner):
    counts = Counter(deficit(row) for row in corner)
    bases = {abs(x) for x in counts if x not in (0, ONE)}
    return len(corner[0]) + counts[ONE] + sum(max(counts[t], counts[-t]) for t in bases)


def extend_permutation(old, rows, cols, size):
    """Return the occupied positions of a permutation extending a scalar corner."""
    occupied_rows = {i for i, j in old}
    occupied_cols = {j for i, j in old}
    unused_rows = [i for i in range(rows) if i not in occupied_rows]
    unused_cols = [j for j in range(cols) if j not in occupied_cols]
    assert size >= rows + cols - len(old)
    return (old + list(zip(unused_rows, range(cols, size)))
            + list(zip(range(rows, size), unused_cols))
            + list(zip(range(rows + len(unused_cols), size),
                       range(cols + len(unused_rows), size))))


def complete(corner):
    """Complete a nonempty rectangular submagic corner at its least size."""
    if not corner or not corner[0] or any(len(row) != len(corner[0]) for row in corner):
        raise ValueError("Supply a nonempty rectangular array")
    columns = list(zip(*corner))
    if not all(submagic_line(row) for row in [*corner, *columns]):
        raise ValueError("Entries must be orthogonal along rows and columns")

    size = minimum_size(corner)
    output = [[0] * size for _ in range(size)]
    fresh_row, fresh_col = len(corner), len(columns)
    blocks = {}
    for i, row in enumerate(corner):
        if not any(row):
            output[i][fresh_col] = ONE
            fresh_col += 1
        for j, entry in enumerate(row):
            if entry == ONE:
                output[i][j] = ONE
            elif entry:
                blocks.setdefault(abs(entry), []).append((i, j, entry))

    for j, column in enumerate(columns):
        if not any(column):
            output[fresh_row][j] = ONE
            fresh_row += 1

    for basis, entries in blocks.items():
        rows = sorted({i for i, j, x in entries})
        cols = sorted({j for i, j, x in entries})
        old_rows, old_cols = len(rows), len(cols)
        row_index = {i: u for u, i in enumerate(rows)}
        col_index = {j: v for v, j in enumerate(cols)}
        patterns = {x: [] for x in (basis, -basis)}
        for i, j, x in entries:
            patterns[x].append((row_index[i], col_index[j]))
        block_size = old_rows + old_cols - min(map(len, patterns.values()))
        rows.extend(range(fresh_row, fresh_row + block_size - old_rows))
        cols.extend(range(fresh_col, fresh_col + block_size - old_cols))
        fresh_row += block_size - old_rows
        fresh_col += block_size - old_cols
        for x, pattern in patterns.items():
            for u, v in extend_permutation(pattern, old_rows, old_cols, block_size):
                i, j = rows[u], cols[v]
                output[i][j] = ONE if output[i][j] else x

    assert fresh_row == fresh_col == size
    return output


@lru_cache(maxsize=None)
def deficits_fit(deficits, width, bases):
    """Independently try all placements of the missing projections in the top border."""
    if width < 0:
        return False
    columns = [[] for _ in range(width)]

    def place(i):
        if i == len(deficits):
            return True
        entry = deficits[i]
        choices = [[(c, entry)] for c in range(width)]
        if entry == ONE:
            for t in range(1, bases + 2):
                choices.extend([(c, t), (d, -t)] for c, d in permutations(range(width), 2))
        for placement in choices:
            if all(all(orthogonal(x, y) for y in columns[c]) for c, x in placement):
                for c, x in placement:
                    columns[c].append(x)
                if place(i + 1):
                    return True
                for c, x in placement:
                    columns[c].pop()
        return False

    # Any extra bases used by zero rows can be identified with one fresh basis:
    # a column that contains a line of such a basis contains no other type.
    return place(0)


def corners(rows, cols, bases):
    alphabet = [0, ONE] + [x for t in range(1, bases + 1) for x in (t, -t)]
    possible_rows = [row for row in product(alphabet, repeat=cols) if submagic_line(row)]

    def extend(old):
        if len(old) == rows:
            yield old
            return
        for row in possible_rows:
            if all(orthogonal(x, y) for prior in old for x, y in zip(row, prior)):
                yield from extend(old + [row])

    yield from extend([])


def check_corners(rows, cols, bases):
    count = 0
    for corner in corners(rows, cols, bases):
        full = complete(corner)
        assert all(tuple(row[:cols]) == old for row, old in zip(full, corner))
        assert all(submagic_line(line) and deficit(line) == 0 for line in [*full, *zip(*full)])
        needs = [deficit(row) for row in corner]
        needs = tuple(sorted((x for x in needs if x != 0), key=str))
        assert not deficits_fit(needs, len(full) - cols - 1, bases)
        count += 1
    return count


# Edge set and orthogonal representation from Cameron, Montanaro, Newman,
# Severini and Winter, Electron. J. Combin. 14 (2007), R81.
EDGES = [
    (1, 2), (1, 3), (1, 11), (1, 12), (1, 16), (2, 3), (2, 4), (2, 13),
    (3, 4), (3, 13), (4, 5), (4, 6), (4, 10), (4, 17), (5, 6), (5, 7),
    (5, 14), (6, 7), (6, 14), (7, 8), (7, 9), (7, 16), (8, 9), (8, 10),
    (8, 13), (9, 10), (9, 13), (10, 11), (10, 12), (10, 17), (11, 12),
    (11, 14), (12, 14), (13, 14), (13, 15), (13, 18), (14, 15), (14, 18),
    (15, 16), (15, 17), (15, 18), (16, 17), (16, 18), (17, 18),
]
VECTORS = [
    (0, 0, 1, -1), (1, 0, 0, 0), (0, 1, 1, 1), (0, 1, 0, -1),
    (0, 0, 1, 0), (1, 1, 0, 1), (1, -1, 0, 0), (0, 0, 0, 1),
    (1, 1, 1, 0), (1, 0, -1, 0), (0, 1, 0, 0), (1, 0, 1, 1),
    (0, 1, -1, 0), (1, 0, 0, -1), (1, 1, 1, 1), (1, 1, -1, -1),
    (1, -1, 1, -1), (1, -1, -1, 1),
]


def coloring(neighbors, colors):
    """Find a proper coloring by backtracking, or return None."""
    assignments = {}

    def search():
        if len(assignments) == len(neighbors):
            return dict(assignments)
        forbidden = {i: {assignments[j] for j in adj if j in assignments}
                     for i, adj in neighbors.items() if i not in assignments}
        i = max(forbidden, key=lambda v: (len(forbidden[v]), len(neighbors[v])))
        for c in range(colors):
            if c in forbidden[i]:
                continue
            assignments[i] = c
            solution = search()
            if solution is not None:
                return solution
            del assignments[i]
        return None

    return search()


def frame(vector):
    a, b, c, d = vector
    return (vector, (-b, a, -d, c), (-c, d, a, -b), (-d, -c, b, a))


def dot(left, right):
    return sum(x * y for x, y in zip(left, right))


def check_graph():
    neighbors = {i: set() for i in range(1, 19)}
    assert len(set(EDGES)) == 44
    for i, j in EDGES:
        neighbors[i].add(j)
        neighbors[j].add(i)
    assert coloring(neighbors, 4) is None
    witness = coloring(neighbors, 5)
    assert witness is not None and all(witness[i] != witness[j] for i, j in EDGES)
    assert set(combinations(range(15, 19), 2)).issubset(EDGES)

    # The four lines at each vertex form an orthogonal basis. Each color
    # assigns orthogonal lines to adjacent vertices. All checks use integers.
    frames = {i: frame(v) for i, v in enumerate(VECTORS, 1)}
    for vectors in frames.values():
        assert all(dot(x, x) > 0 for x in vectors)
        assert all(dot(x, y) == 0 for x, y in combinations(vectors, 2))
    for i, j in EDGES:
        assert all(dot(x, y) == 0 for x, y in zip(frames[i], frames[j]))


if __name__ == "__main__":
    for rows, cols in [(2, 3), (3, 2), (3, 3)]:
        count = check_corners(rows, cols, 3)
        print(f"{rows} x {cols}, three complementary pairs: {count} corners checked")
    check_graph()
    print("18 vertices, 44 edges: completion sizes 23 over M_18 and 22 after tensoring with I_4")
