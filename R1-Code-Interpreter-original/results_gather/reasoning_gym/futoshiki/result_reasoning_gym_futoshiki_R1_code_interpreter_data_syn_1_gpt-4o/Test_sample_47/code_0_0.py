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
        [4, 0, 0, 2, 0],
        [0, 0, 5, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [2, 0, 3, 0, 0]
    ]
    
    # Constraints: (row1, col1, row2, col2, operator)
    constraints = [
        (0, 1, 0, 2, '>'),
        (1, 0, 1, 1, '>'),
        (1, 1, 1, 2, '<'),
        (1, 2, 1, 3, '>'),
        (2, 3, 2, 4, '>'),
        (3, 1, 3, 2, '<'),
        (3, 3, 3, 4, '<'),
        (0, 1, 1, 1, '>'),
        (1, 1, 2, 1, '<'),
        (1, 3, 2, 3, '>'),
        (2, 2, 3, 2, '<')
    ]
    
    # Try all permutations for each row
    for perm1 in permutations(range(1, 6)):
        if perm1[0] == 4 and perm1[3] == 2:
            grid[0] = list(perm1)
            for perm2 in permutations(range(1, 6)):
                if perm2[2] == 5:
                    grid[1] = list(perm2)
                    for perm3 in permutations(range(1, 6)):
                        grid[2] = list(perm3)
                        for perm4 in permutations(range(1, 6)):
                            grid[3] = list(perm4)
                            for perm5 in permutations(range(1, 6)):
                                if perm5[0] == 2 and perm5[2] == 3:
                                    grid[4] = list(perm5)
                                    if is_valid(grid, constraints):
                                        for row in grid:
                                            print(row)
                                        return

solve_futoshiki()