"""Exact checks for the fixed-rank diagonal algebras in Section 6.

In a symbolic array, t and -t denote the free projection numbered t and
its complement. Zero and ONE denote scalar entries.
"""

from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations, product

from diagonals import difference, multiply, projector, transpose


ONE = "I"


def universal_array(n, rank, marks):
    """Give a symbolic realization for every free set of marked positions."""
    if not (1 <= marks <= n and marks <= rank <= min(n, 2 * (n - marks))):
        raise ValueError("The marked projections are not universal at these parameters")
    u = [[0] * n for _ in range(n)]
    blocks = min(marks, rank - marks)
    for i in range(marks):
        u[i][i] = i + 1
    for i in range(blocks):
        j = marks + i
        u[i][j] = u[j][i] = -(i + 1)
        u[j][j] = i + 1
    ports = [(j, is_column) for j in range(marks + blocks, n) for is_column in (True, False)]
    for i, (j, is_column) in zip(range(blocks, marks), ports):
        row, col = (i, j) if is_column else (j, i)
        u[row][col] = -(i + 1)
    for j in range(2 * marks, rank):
        u[j][j] = ONE
    return u


def check_symbolic(u, rank, marks):
    for line in [*u, *zip(*u)]:
        used = [x for x in line if x != 0]
        if len(used) > 1:
            assert len(used) == 2 and ONE not in used and used[0] == -used[1]
    entries = Counter(x for row in u for x in row)
    assert all(u[i][i] == i + 1 for i in range(marks))
    assert all(entries[t] == entries[-t] for t in range(1, marks + 1))
    assert entries[ONE] + sum(entries[-t] for t in range(1, marks + 1)) == rank


def check_threshold(max_n):
    realized = obstructions = 0
    for n in range(1, max_n + 1):
        for marks in range(1, n + 1):
            # The deficit lines have vectors (1, slope); the vector state
            # at (1, 0) sees slope^2/(1+slope^2) on their complements.
            slopes = [Fraction(t, 4 * marks ** 2) for t in range(1, marks + 1)]
            residual = sum(s ** 2 / (1 + s ** 2) for s in slopes)
            assert len(set(slopes)) == marks and residual < 1
            for rank in range(n + 1):
                if marks > rank:
                    obstructions += 1  # All marked projections equal to one.
                elif rank > 2 * (n - marks):
                    assert 2 * (n - marks) + residual < rank
                    obstructions += 1
                else:
                    check_symbolic(universal_array(n, rank, marks), rank, marks)
                    realized += 1
    return realized, obstructions


def check_patterns(n):
    # Keep fixed-point sets rather than their cardinalities.
    for rank in range(n + 1):
        actual = set()
        for domain in combinations(range(n), rank):
            for image in permutations(range(n), rank):
                actual.add(frozenset(a for a, i in zip(domain, image) if a == i))
        sizes = range(rank + 1) if rank < n else [k for k in range(n + 1) if k != n - 1]
        expected = {frozenset(s) for k in sizes for s in combinations(range(n), k)}
        assert actual == expected


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def check_three_point_fibers():
    # Sample the diagonal sign-change orbits in the Pauli model, including
    # fibers with repeated lines and fibers with four orthogonal lines.
    signs = [(1, 1, 1, 1), (1, 1, -1, -1), (1, -1, 1, -1), (1, -1, -1, 1)]
    for coords in product(range(3), repeat=4):
        if not any(coords):
            continue
        lines = [tuple(a * b for a, b in zip(s, coords)) for s in signs]
        for a, b, c, d in permutations(lines):
            assert abs(dot(a, b)) == abs(dot(c, d))
        # The component of the killed line is excluded from any reducing
        # summand killing that projection. A surviving component has at most two labels.
        pending = set(range(4))
        while pending:
            reached = {min(pending)}
            while True:
                adjoining = {i for i in pending if any(dot(lines[i], lines[j]) for j in reached)}
                if adjoining <= reached:
                    break
                reached |= adjoining
            pending -= reached
            assert 3 in reached or len(reached) <= 2


def check_rank_two_examples():
    # Bloch vectors of a regular tetrahedron, in units of 1/sqrt(3).
    bloch = [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]
    assert all(sum(col) == 0 for col in zip(*bloch))
    assert all(dot(v, v) == 3 for v in bloch)
    for x, y in combinations(bloch, 2):
        overlap = Fraction(1, 2) * (1 + Fraction(dot(x, y), 3))
        assert overlap == Fraction(1, 3)

    # Rational instance of the four-point array with q and r orthogonal.
    unit = [[int(i == j) for j in range(3)] for i in range(3)]
    zero = [[0] * 3 for _ in range(3)]
    p, q, r = [projector(v) for v in ([1, 2, 3], [1, 0, 0], [0, 1, 0])]
    u = [[p, zero, zero, difference(unit, p)], [zero, q, zero, zero],
         [zero, zero, r, zero], [zero, difference(difference(unit, q), r), zero, zero]]
    for line in [*u, *zip(*u)]:
        for a in line:
            assert multiply(a, a) == a and list(map(list, transpose(a))) == a
        assert all(multiply(a, b) == zero for a, b in combinations(line, 2))
    for i, j in product(range(3), repeat=2):
        assert sum(a[i][j] for row in u for a in row) == 2 * unit[i][j]

    pairs = {frozenset(p) for p in combinations(range(3), 2)}
    for length in range(1, 5):
        covers = []
        for word in product(range(3), repeat=length):
            edges = list(zip(word, word[1:]))
            if all(a != b for a, b in edges) and {frozenset(e) for e in edges} == pairs:
                covers.append(word)
        assert len(covers) == (6 if length == 4 else 0)


def check_four_letter_witnesses():
    # For the cyclic triple of tetrahedral projections c=(1+2w)/9;
    # conjugating w to -1-w changes its sign.
    # <v_1,v_2><v_2,v_3><v_3,v_1> gives (1+2w)/9.
    c = (Fraction(1, 9), Fraction(2, 9))
    conjugate = (c[0] - c[1], -c[1])
    assert c != (0, 0) and conjugate == (-c[0], -c[1])
    tetra_gram = [[Fraction(1) if i == j else Fraction(1, 3) for j in range(3)]
                  for i in range(3)]
    signs = [[0, 1, 1, 1], [-1, 0, 1, -1], [-1, -1, 0, 1], [-1, 1, -1, 0]]
    assert all(dot(x, y) == (3 if i == j else 0) for i, x in enumerate(signs)
               for j, y in enumerate(signs))
    # f_i and g_i are orthonormal bases of their respective spans.
    # Their mixed inner products are squares of entries of signs/sqrt(3).
    def inner(i, a, j, b):
        if a == b:
            return int(i == j)
        return Fraction(signs[i][j] ** 2, 3)

    def coefficient(labels, at):
        # Compute <f_at, product(p_label) g_at> by expanding each rank-two projection.
        result = 0
        for types in product((0, 1), repeat=len(labels)):
            sites = [(at, 0), *zip(labels, types), (at, 1)]
            term = 1
            for (i, a), (j, b) in zip(sites, sites[1:]):
                term *= inner(i, a, j, b)
            result += term
        return result

    second_gram = []
    for at in range(3):
        row = []
        for i in range(3):
            j, k = (i + 1) % 3, (i + 2) % 3
            entry = 27 * coefficient([i, j, k, i], at)
            assert entry == 27 * coefficient([i, k, j, i], at)
            assert entry == (1 if at == i else Fraction(1, 9))
            row.append(entry)
        second_gram.append(row)
    for gram in (tetra_gram, second_gram):
        # Strict diagonal dominance certifies invertibility in exact arithmetic.
        assert all(row[i] > sum(abs(x) for j, x in enumerate(row) if j != i)
                   for i, row in enumerate(gram))


if __name__ == "__main__":
    for order in range(1, 8):
        check_patterns(order)
    models, obstructions = check_threshold(20)
    check_three_point_fibers()
    check_rank_two_examples()
    check_four_letter_witnesses()
    p, q, r = (projector(v) for v in ([1, 0], [1, 1], [0, 1]))
    assert multiply(multiply(p, q), r) == [[0, Fraction(1, 2)], [0, 0]]
    print("All classical fixed-point sets through seven points checked.")
    print(f"Through N=20: {models} symbolic rank constructions and {obstructions} obstruction checks.")
    print("Three-point Pauli fibers, four-point witnesses and reduced words through degree four checked.")
