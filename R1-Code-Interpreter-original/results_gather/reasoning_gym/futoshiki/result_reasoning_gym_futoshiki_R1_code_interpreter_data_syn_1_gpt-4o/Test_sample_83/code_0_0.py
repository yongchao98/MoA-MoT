from itertools import permutations

# Initial grid with constraints
grid = [
    [0, 3, 0, 6, 5, 0],
    [0, 0, 0, 2, 0, 6],
    [0, 0, 4, 0, 6, 0],
    [0, 0, 1, 0, 0, 0],
    [0, 0, 0, 0, 0, 5],
    [0, 5, 6, 3, 0, 0]
]

# Inequality constraints
horizontal_constraints = [
    [None, None, None, None, None],
    [None, None, None, None, None],
    [None, None, None, None, '>'],
    [None, None, None, None, None],
    ['>', None, None, None, None],
    [None, None, None, None, None]
]

vertical_constraints = [
    [None, None, None, None, None, None],
    [None, None, None, None, None, None],
    [None, None, None, None, None, None],
    [None, None, None, None, None, None],
    [None, None, None, None, None, None],
    [None, None, None, None, None, None]
]

# Add vertical constraints
vertical_constraints[2][5] = '>'
vertical_constraints[4][1] = '<'

def is_valid(grid, row, col, num):
    # Check row and column uniqueness
    for i in range(6):
        if grid[row][i] == num or grid[i][col] == num:
            return False
    
    # Check horizontal constraints
    if col > 0 and horizontal_constraints[row][col - 1]:
        if horizontal_constraints[row][col - 1] == '>' and not (grid[row][col - 1] > num):
            return False
        if horizontal_constraints[row][col - 1] == '<' and not (grid[row][col - 1] < num):
            return False
    if col < 5 and horizontal_constraints[row][col]:
        if horizontal_constraints[row][col] == '>' and not (num > grid[row][col + 1]):
            return False
        if horizontal_constraints[row][col] == '<' and not (num < grid[row][col + 1]):
            return False
    
    # Check vertical constraints
    if row > 0 and vertical_constraints[row - 1][col]:
        if vertical_constraints[row - 1][col] == '>' and not (grid[row - 1][col] > num):
            return False
        if vertical_constraints[row - 1][col] == '<' and not (grid[row - 1][col] < num):
            return False
    if row < 5 and vertical_constraints[row][col]:
        if vertical_constraints[row][col] == '>' and not (num > grid[row + 1][col]):
            return False
        if vertical_constraints[row][col] == '<' and not (num < grid[row + 1][col]):
            return False
    
    return True

def solve(grid):
    for row in range(6):
        for col in range(6):
            if grid[row][col] == 0:
                for num in range(1, 7):
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