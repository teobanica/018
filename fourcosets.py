"""Check four-point differences and the four-layer examples by exact arithmetic."""

from itertools import combinations, product

from quotients import analyze, band, quotient


def check_group(moduli):
    """Check every four-point set up to translation to a set containing zero."""
    group = list(product(*(range(n) for n in moduli)))
    index = {g: i for i, g in enumerate(group)}
    size = len(group)
    add = [[index[tuple((a + b) % n for a, b, n in zip(g, h, moduli))]
            for h in group] for g in group]
    neg = [index[tuple(-a % n for a, n in zip(g, moduli))] for g in group]
    forbidden = {x for x in range(1, size)
                 if add[x][x] == 0 or add[add[x][x]][x] == 0}
    progressions = set()
    for start, step in product(range(size), repeat=2):
        points = [start]
        for _ in range(3):
            points.append(add[points[-1]][step])
        if len(set(points)) == 4 and 0 in points:
            progressions.add(frozenset(points))

    count = 0
    for others in combinations(range(1, size), 3):
        points = (0, *others)
        diffs = {add[a][neg[b]] for a, b in product(points, repeat=2)}
        if len(diffs) != 7 or diffs & forbidden:
            continue
        if any(all(add[d][t] in diffs for d in diffs) for t in diffs if t != 0):
            continue
        assert frozenset(points) in progressions, (moduli, points)
        count += 1
    return count


def check_orders():
    """Read the stabilizer and step order directly from the matrix entries."""
    for order in (14, 21, 28, 35, 42):
        rows, level = band(4, 2, phase_order=2 * order)
        assert analyze(rows, level)["commuting"]
        diffs = {quotient(a, b, level) for a, b in product(rows, repeat=2)}
        stabilizer = {h for h in diffs
                      if all(quotient(d, h, level) in diffs for d in diffs)}
        assert len(diffs) == 14 and stabilizer == {tuple(rows[0]), tuple(rows[1])}
        step = rows[2]
        first_return = next(t for t in range(1, level + 1)
                            if tuple(t * a % level for a in step) in stabilizer)
        assert first_return == order


if __name__ == "__main__":
    shapes = [(n,) for n in range(4, 36)] + [(2, 4), (3, 3), (4, 4), (2, 14),
                                               (3, 7), (5, 5), (2, 2, 7)]
    count = sum(check_group(shape) for shape in shapes)
    check_orders()
    print(f"{count} four-point sets passing the difference tests in {len(shapes)} groups.")
    print("Seven-column fibers and step orders 14, 21, 28, 35 and 42 verified.")
