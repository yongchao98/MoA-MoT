from itertools import permutations

# Initial grid setup with given numbers and constraints
grid = [
    [0, 0, 4, 0, 0, 0, 5, 3],
    [0, 5, 3, 4, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 5, 7, 0],
    [0, 0, 0, 0, 0, 2, 0, 0],
    [0, 0, 2, 5, 0, 0, 0, 6],
    [1, 0, 8, 0, 0, 7, 4, 0],
    [0, 0, 0, 0, 4, 0, 0, 0],
    [0, 0, 5, 6, 3, 0, 0, 8]
]

# Constraints: (row, col, direction, row2, col2)
# direction: 1 for '<', -1 for '>', 2 for '∧', -2 for '∨'
constraints = [
    (0, 2, -1, 0, 3), (0, 6, 1, 0, 7),
    (1, 0, 2, 0, 0), (1, 2, 2, 0, 2), (1, 7, 1, 1, 6),
    (2, 5, 2, 1, 5), (2, 6, -1, 2, 7),
    (3, 2, 1, 3, 1), (3, 5, 2, 2, 5), (3, 7, 2, 2, 7),
    (4, 0, 2, 3, 0), (4, 3, 2, 3, 3), (4, 5, 2, 3, 5),
    (5, 1, 1, 5, 0), (5, 3, 1, 5, 2), (5, 5, 2, 4, 5),
    (6, 0, -1, 6, 1), (6, 2, 1, 6, 3), (6, 4, -1, 6, 5), (6, 6, 1, 6, 5),
    (7, 3, -1, 7, 4), (7, 4, 1, 7, 5)
]

def is_valid(grid, constraints):
    # Check all constraints
    for r1, c1, direction, r2, c2 in constraints:
        if direction == 1 and not (grid[r1][c1] < grid[r2][c2]):
            return False
        if direction == -1 and not (grid[r1][c1] > grid[r2][c2]):
            return False
        if direction == 2 and not (grid[r1][c1] < grid[r2][c2]):
            return False
        if direction == -2 and not (grid[r1][c1] > grid[r2][c2]):
            return False
    return True

def solve(grid, constraints):
    # Find the first empty cell
    for r in range(8):
        for c in range(8):
            if grid[r][c] == 0:
                for num in range(1, 9):
                    grid[r][c] = num
                    if is_valid(grid, constraints) and solve(grid, constraints):
                        return True
                    grid[r][c] = 0
                return False
    return True

solve(grid, constraints)

# Print the solved grid
for row in grid:
    print(row)