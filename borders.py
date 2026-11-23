"""Check rectangular magic completions and transpose identities over the rationals."""

from itertools import combinations, permutations, product

from diagonals import difference, multiply, projector, transpose, verify_magic, word


def add_matrices(matrices, dimension):
    return [[sum(matrix[i][j] for matrix in matrices) for j in range(dimension)]
            for i in range(dimension)]


def magic_border(entries):
    """Complete an m by n array of real orthogonal projections to size m+n."""
    m, n = len(entries), len(entries[0])
    dimension = len(entries[0][0])
    identity = [[int(a == b) for b in range(dimension)] for a in range(dimension)]
    zero = add_matrices([], dimension)
    row_deficits = [difference(identity, add_matrices(row, dimension)) for row in entries]
    col_deficits = [difference(identity, add_matrices(col, dimension)) for col in transpose(entries)]
    upper = [list(row) + [row_deficits[i] if i == j else zero for j in range(m)]
             for i, row in enumerate(entries)]
    lower = [[col_deficits[i] if i == j else zero for j in range(n)] + list(col)
             for i, col in enumerate(transpose(entries))]
    return upper + lower


def check_scalar_corners():
    checked = 0
    for m, n in product(range(1, 4), repeat=2):
        for rank in range(min(m, n) + 1):
            for source in combinations(range(n), rank):
                for target in permutations(range(m), rank):
                    arrows = set(zip(target, source))
                    entries = [[[[int((i, j) in arrows)]] for j in range(n)] for i in range(m)]
                    full = magic_border(entries)
                    verify_magic(full)
                    assert [row[:n] for row in full[:m]] == entries
                    checked += 1
        for size in range(max(m, n), m + n + 1):
            ranks = {sum(i < m for i in perm[:n]) for perm in permutations(range(size))}
            assert ranks == set(range(max(0, m + n - size), min(m, n) + 1))
    return checked


def inverse_expression(entries, indices):
    """Evaluate m3(S x id x S) Delta2 on a word, reversing the outer products."""
    size, dimension = len(entries), len(entries[0][0])
    terms = []
    for pairs in product(product(range(size), repeat=2), repeat=len(indices)):
        left = [entries[k][i] for (i, j), (k, l) in zip(indices, pairs)]
        middle = [entries[k][l] for k, l in pairs]
        right = [entries[j][l] for (i, j), (k, l) in zip(indices, pairs)]
        terms.append(word(*reversed(left), *middle, *reversed(right)))
    return add_matrices(terms, dimension)


def check_operator_corners():
    p, q = projector([1, 0]), projector([1, 1])
    zero = [[0, 0], [0, 0]]
    verify_magic(magic_border([[p, zero, zero], [zero, q, zero]]))
    entries = [[p, zero], [zero, q]]
    for i, j in product(range(2), repeat=2):
        assert inverse_expression(entries, [(i, j)]) == entries[j][i]
    left = inverse_expression(entries, [(0, 0), (1, 1)])
    right = multiply(q, p)
    assert left != right
    assert add_matrices([left, left], 2) == right
    assert word(p, q, p, q) != multiply(p, q)


if __name__ == "__main__":
    count = check_scalar_corners()
    check_operator_corners()
    print(f"Verified {count} scalar completions, all attainable corner ranks through 3 by 3,")
    print("a noncommuting rectangular completion, and the degree-two transpose counterexample.")
