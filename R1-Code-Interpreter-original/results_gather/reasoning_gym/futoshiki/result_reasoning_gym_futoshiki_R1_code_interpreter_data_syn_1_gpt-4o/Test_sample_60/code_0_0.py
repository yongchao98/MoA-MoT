from itertools import permutations

def is_valid(grid, constraints):
    for (r1, c1, r2, c2, op) in constraints:
        if op == '<' and not (grid[r1][c1] < grid[r2][c2]):
            return False
        if op == '>' and not (grid[r1][c1] > grid[r2][c2]):
            return False
        if op == '∧' and not (grid[r1][c1] < grid[r2][c2]):
            return False
        if op == '∨' and not (grid[r1][c1] > grid[r2][c2]):
            return False
    return True

def solve_futoshiki():
    # Initial grid with known values
    grid = [
        [0, 0, 5, 0, 0, 0, 0],
        [0, 0, 0, 1, 0, 0, 0],
        [0, 1, 0, 0, 0, 0, 2],
        [0, 0, 1, 0, 2, 0, 4],
        [0, 0, 0, 0, 0, 0, 1],
        [0, 0, 0, 0, 5, 0, 0],
        [0, 0, 0, 0, 7, 0, 0]
    ]

    # Constraints (r1, c1, r2, c2, op)
    constraints = [
        (0, 5, 0, 6, '<'),
        (2, 2, 2, 3, '<'),
        (3, 2, 3, 3, '<'),
        (5, 3, 5, 4, '>'),
        (6, 0, 6, 1, '<'),
        (6, 5, 6, 6, '<'),
        (1, 0, 2, 0, '∧'),
        (2, 2, 3, 2, '∨')
    ]

    # Try permutations for each row
    for perm1 in permutations(range(1, 8)):
        if perm1[2] != 5:
            continue
        grid[0] = list(perm1)
        for perm2 in permutations(range(1, 8)):
            if perm2[3] != 1:
                continue
            grid[1] = list(perm2)
            for perm3 in permutations(range(1, 8)):
                if perm3[1] != 1 or perm3[6] != 2:
                    continue
                grid[2] = list(perm3)
                for perm4 in permutations(range(1, 8)):
                    if perm4[2] != 1 or perm4[4] != 2 or perm4[6] != 4:
                        continue
                    grid[3] = list(perm4)
                    for perm5 in permutations(range(1, 8)):
                        if perm5[6] != 1:
                            continue
                        grid[4] = list(perm5)
                        for perm6 in permutations(range(1, 8)):
                            if perm6[4] != 5:
                                continue
                            grid[5] = list(perm6)
                            for perm7 in permutations(range(1, 8)):
                                if perm7[4] != 7:
                                    continue
                                grid[6] = list(perm7)
                                if is_valid(grid, constraints):
                                    return grid

solution = solve_futoshiki()
for row in solution:
    print(row)