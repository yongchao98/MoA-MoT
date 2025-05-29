from itertools import permutations

def is_valid(grid, constraints):
    for (r1, c1, r2, c2, op) in constraints:
        if op == '>' and not (grid[r1][c1] > grid[r2][c2]):
            return False
        if op == '<' and not (grid[r1][c1] < grid[r2][c2]):
            return False
    return True

def solve_futoshiki():
    # Initial grid with known values
    grid = [
        [0, 0, 0, 0, 3, 4],
        [0, 0, 5, 0, 0, 0],
        [5, 2, 0, 0, 0, 3],
        [0, 0, 0, 2, 0, 0],
        [0, 4, 0, 0, 0, 2],
        [0, 0, 0, 0, 2, 0]
    ]

    # Constraints: (row1, col1, row2, col2, operator)
    constraints = [
        (0, 2, 0, 3, '>'),
        (1, 1, 1, 2, '>'),
        (2, 1, 2, 2, '>'),
        (2, 4, 2, 5, '>'),
        (3, 0, 3, 1, '<'),
        (3, 3, 3, 4, '<'),
        (3, 4, 3, 5, '<'),
        (4, 5, 4, 4, '>'),
        (5, 0, 5, 1, '>'),
        (5, 1, 5, 2, '<'),
        (5, 4, 5, 5, '<'),
        (0, 0, 1, 0, 'v'),
        (1, 2, 2, 2, 'v'),
        (2, 5, 3, 5, 'v'),
        (3, 1, 4, 1, 'v'),
        (4, 5, 5, 5, 'v')
    ]

    # Try all permutations of numbers 1 to 6 for each row
    for perm in permutations(range(1, 7)):
        grid[0][:3] = perm[:3]
        grid[1][:2] = perm[3:5]
        grid[2][2] = perm[5]
        if is_valid(grid, constraints):
            for perm2 in permutations(range(1, 7)):
                grid[3][:3] = perm2[:3]
                grid[4][:2] = perm2[3:5]
                grid[5][2] = perm2[5]
                if is_valid(grid, constraints):
                    for perm3 in permutations(range(1, 7)):
                        grid[0][3:5] = perm3[:2]
                        grid[1][2:4] = perm3[2:4]
                        grid[2][3:5] = perm3[4:]
                        if is_valid(grid, constraints):
                            for perm4 in permutations(range(1, 7)):
                                grid[3][3:5] = perm4[:2]
                                grid[4][2:4] = perm4[2:4]
                                grid[5][3:5] = perm4[4:]
                                if is_valid(grid, constraints):
                                    for perm5 in permutations(range(1, 7)):
                                        grid[0][5] = perm5[0]
                                        grid[1][5] = perm5[1]
                                        grid[2][5] = perm5[2]
                                        grid[3][5] = perm5[3]
                                        grid[4][5] = perm5[4]
                                        grid[5][5] = perm5[5]
                                        if is_valid(grid, constraints):
                                            for perm6 in permutations(range(1, 7)):
                                                grid[0][0] = perm6[0]
                                                grid[1][0] = perm6[1]
                                                grid[2][0] = perm6[2]
                                                grid[3][0] = perm6[3]
                                                grid[4][0] = perm6[4]
                                                grid[5][0] = perm6[5]
                                                if is_valid(grid, constraints):
                                                    for perm7 in permutations(range(1, 7)):
                                                        grid[0][1] = perm7[0]
                                                        grid[1][1] = perm7[1]
                                                        grid[2][1] = perm7[2]
                                                        grid[3][1] = perm7[3]
                                                        grid[4][1] = perm7[4]
                                                        grid[5][1] = perm7[5]
                                                        if is_valid(grid, constraints):
                                                            for perm8 in permutations(range(1, 7)):
                                                                grid[0][2] = perm8[0]
                                                                grid[1][2] = perm8[1]
                                                                grid[2][2] = perm8[2]
                                                                grid[3][2] = perm8[3]
                                                                grid[4][2] = perm8[4]
                                                                grid[5][2] = perm8[5]
                                                                if is_valid(grid, constraints):
                                                                    return grid

solution = solve_futoshiki()
for row in solution:
    print(row)