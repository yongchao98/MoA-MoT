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
        if perm1[3] != 3:  # Respect the initial value
            continue
        grid[0] = list(perm1)
        if not is_valid(grid, constraints):
            continue
        for perm2 in permutations(range(1, 6)):
            grid[1] = list(perm2)
            if not is_valid(grid, constraints):
                continue
            for perm3 in permutations(range(1, 6)):
                if perm3[0] != 2:  # Respect the initial value
                    continue
                grid[2] = list(perm3)
                if not is_valid(grid, constraints):
                    continue
                for perm4 in permutations(range(1, 6)):
                    if perm4[1] != 1 or perm4[3] != 5 or perm4[4] != 4:  # Respect the initial values
                        continue
                    grid[3] = list(perm4)
                    if not is_valid(grid, constraints):
                        continue
                    for perm5 in permutations(range(1, 6)):
                        if perm5[0] != 1:  # Respect the initial value
                            continue
                        grid[4] = list(perm5)
                        if not is_valid(grid, constraints):
                            continue
                        # Check columns
                        if all(len(set(col)) == 5 for col in zip(*grid)):
                            # Check if the grid is valid with all constraints
                            if is_valid(grid, constraints):
                                return grid

solution = solve_futoshiki()
if solution:
    for row in solution:
        print(row)
else:
    print("No solution found")