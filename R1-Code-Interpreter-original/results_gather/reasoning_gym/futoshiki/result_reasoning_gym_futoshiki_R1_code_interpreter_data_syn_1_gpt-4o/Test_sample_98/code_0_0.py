from itertools import permutations

# Initial grid with known values and constraints
grid = [
    [0, 0, 1, 0, 0, 0, 0, 5],
    [0, 0, 0, 0, 2, 0, 0, 0],
    [0, 0, 2, 0, 8, 7, 0, 0],
    [8, 0, 0, 0, 0, 4, 0, 1],
    [0, 0, 3, 4, 0, 1, 0, 7],
    [4, 0, 0, 7, 0, 0, 0, 0],
    [5, 0, 0, 0, 0, 2, 6, 0],
    [0, 0, 4, 5, 7, 0, 0, 0]
]

# Constraints
horizontal_constraints = [
    [None, None, None, None, '<', None, None],
    [None, None, '>', None, None, '<', None],
    ['>', '>', None, None, None, '>', None],
    [None, '>', None, None, '>', '>', None],
    [None, None, None, '>', None, None, '<'],
    ['>', '<', None, None, '<', None, '<'],
    [None, None, None, None, None, None, None]
]

vertical_constraints = [
    [None, None, None, None, None, None, None, None],
    ['∨', None, None, None, '∨', None, '∨', None],
    [None, None, None, None, None, None, None, None],
    ['∧', '∧', '∧', None, None, '∨', None, None],
    [None, '∨', None, None, None, None, '∧', None],
    ['∧', None, None, None, None, None, None, None],
    [None, '∨', None, None, '∧', '∧', None, None]
]

def is_valid(grid, row, col, num):
    # Check row and column uniqueness
    for i in range(8):
        if grid[row][i] == num or grid[i][col] == num:
            return False

    # Check horizontal constraints
    if col > 0 and horizontal_constraints[row][col - 1]:
        if horizontal_constraints[row][col - 1] == '<' and not (grid[row][col - 1] < num):
            return False
        if horizontal_constraints[row][col - 1] == '>' and not (grid[row][col - 1] > num):
            return False
    if col < 7 and horizontal_constraints[row][col]:
        if horizontal_constraints[row][col] == '<' and not (num < grid[row][col + 1]):
            return False
        if horizontal_constraints[row][col] == '>' and not (num > grid[row][col + 1]):
            return False

    # Check vertical constraints
    if row > 0 and vertical_constraints[row - 1][col]:
        if vertical_constraints[row - 1][col] == '∧' and not (grid[row - 1][col] > num):
            return False
        if vertical_constraints[row - 1][col] == '∨' and not (grid[row - 1][col] < num):
            return False
    if row < 7 and vertical_constraints[row][col]:
        if vertical_constraints[row][col] == '∧' and not (num > grid[row + 1][col]):
            return False
        if vertical_constraints[row][col] == '∨' and not (num < grid[row + 1][col]):
            return False

    return True

def solve(grid):
    for row in range(8):
        for col in range(8):
            if grid[row][col] == 0:
                for num in range(1, 9):
                    if is_valid(grid, row, col, num):
                        grid[row][col] = num
                        if solve(grid):
                            return True
                        grid[row][col] = 0
                return False
    return True

solve(grid)

# Print the solved grid
for row in grid:
    print(row)