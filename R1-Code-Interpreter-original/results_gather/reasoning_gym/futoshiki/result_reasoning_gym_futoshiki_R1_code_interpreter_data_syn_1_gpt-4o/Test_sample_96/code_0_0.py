from itertools import permutations

def is_valid(grid, constraints):
    for (r1, c1, r2, c2, op) in constraints:
        if op == '<' and not (grid[r1][c1] < grid[r2][c2]):
            return False
        if op == '>' and not (grid[r1][c1] > grid[r2][c2]):
            return False
    return True

def solve_futoshiki():
    # Initial grid with known values
    grid = [
        [0, 0, 0, 3, 0],
        [0, 0, 0, 0, 0],
        [2, 0, 0, 0, 0],
        [0, 1, 0, 5, 4],
        [1, 0, 0, 0, 0]
    ]

    # Constraints: (row1, col1, row2, col2, operator)
    constraints = [
        (0, 3, 0, 4, '<'),
        (1, 0, 1, 1, '<'),
        (2, 1, 2, 2, '<'),
        (4, 3, 4, 4, '>')
    ]

    # Try all permutations of numbers 1 to 5 for each row
    for perm1 in permutations(range(1, 6)):
        grid[0] = list(perm1)
        if not is_valid(grid, constraints):
            continue
        for perm2 in permutations(range(1, 6)):
            grid[1] = list(perm2)
            if not is_valid(grid, constraints):
                continue
            for perm3 in permutations(range(1, 6)):
                grid[2] = list(perm3)
                if not is_valid(grid, constraints):
                    continue
                for perm5 in permutations(range(1, 6)):
                    grid[4] = list(perm5)
                    if not is_valid(grid, constraints):
                        continue
                    # Check columns
                    if all(len(set(col)) == 5 for col in zip(*grid)):
                        # Check if the grid is valid with all constraints
                        if is_valid(grid, constraints):
                            return grid

solution = solve_futoshiki()
for row in solution:
    print(row)