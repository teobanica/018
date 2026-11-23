"""Exact rational checks of the diagonal relation in four and five points."""

from fractions import Fraction
from itertools import permutations


def transpose(a):
    return list(zip(*a))


def multiply(a, b):
    return [[sum(x * y for x, y in zip(row, col)) for col in transpose(b)] for row in a]


def word(*factors):
    result = factors[0]
    for factor in factors[1:]:
        result = multiply(result, factor)
    return result


def difference(a, b):
    return [[x - y for x, y in zip(row, other)] for row, other in zip(a, b)]


def projector(v):
    norm2 = sum(x * x for x in v)
    return [[Fraction(x * y, norm2) for y in v] for x in v]


def witness(p, q, r):
    return difference(word(p, q, p, r, p), word(p, r, p, q, p))


def verify_magic(u):
    n, d = len(u), len(u[0][0])
    for entries in (*u, *transpose(u)):
        for i in range(d):
            for j in range(d):
                assert sum(a[i][j] for a in entries) == (i == j)
        for a in entries:
            assert a == multiply(a, a)
            assert list(map(list, transpose(a))) == a


def order_five():
    dim = 3
    identity = [[int(i == j) for j in range(dim)] for i in range(dim)]
    zero = [[0] * dim for _ in range(dim)]
    complement = projector([1, 1, 1])
    p = difference(identity, complement)
    u = [[zero for _ in range(5)] for _ in range(5)]
    u[0][0] = u[1][1] = p
    u[0][1] = u[1][0] = complement
    perms = [(0, 2, 1), (2, 1, 0), (1, 0, 2)]
    for i in range(dim):
        for j in range(dim):
            u[2 + i][2 + j] = [[int(h == l and perms[h][j] == i) for l in range(dim)]
                                for h in range(dim)]
    verify_magic(u)
    result = witness(u[0][0], u[2][2], u[3][3])
    scaled = [[9 * x for x in row] for row in result]
    assert scaled == [[0, -1, 1], [1, 0, -1], [-1, 1, 0]]
    return [[int(x) for x in row] for row in scaled]


def order_four():
    # Pauli diagonal lines are quaternion conjugates, hence coordinate sign changes.
    signs = [(1, 1, 1, 1), (1, 1, -1, -1), (1, -1, 1, -1), (1, -1, -1, 1)]
    for coords in ([1, 2, 3, 4], [0, 1, 1, 3], [4, 0, 0, 0]):
        points = [projector([a * b for a, b in zip(sign, coords)]) for sign in signs]
        for p, q, r in permutations(points, 3):
            assert all(x == 0 for row in witness(p, q, r) for x in row)


def fixed_patterns(n):
    patterns = {tuple(int(perm[i] != i) for i in range(n)) for perm in permutations(range(n))}
    assert len(patterns) == 2 ** n - n
    assert all(sum(y) != 1 for y in patterns)
    return len(patterns)


if __name__ == "__main__":
    order_four()
    print("Pauli diagonal checks passed")
    print("Nine times the order-five witness:")
    for row in order_five():
        print(*row)
    for n in range(1, 9):
        print(n, "points,", fixed_patterns(n), "classical diagonal characters")
