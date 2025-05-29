# Initial grid with constraints
grid = [
    [0, 0, 0, 1, 0, 0, 4],
    [2, 0, 0, 0, 0, 0, 5],
    [7, 0, 3, 0, 0, 0, 2],
    [0, 0, 0, 4, 0, 5, 0],
    [0, 0, 0, 0, 5, 0, 0],
    [0, 7, 4, 0, 6, 0, 0],
    [0, 3, 0, 5, 2, 0, 7]
]

# Inequality constraints
horizontal_constraints = [
    [None, None, '<', None, '<', None],
    [None, '>', '<', None, None, None],
    [None, '>', None, None, None, None],
    [None, '>', None, '>', '>', None],
    [None, None, None, None, None, None],
    [None, None, None, '>', None, None],
    [None, None, None, None, None, None]
]

vertical_constraints = [
    [None, None, None, None, None, None, None],
    [None, None, None, None, None, None, None],
    [None, None, None, None, None, None, None],
    ['∨', None, None, None, None, None, None],
    [None, None, None, None, None, None, None],
    [None, None, None, None, None, None, None]
]

def is_valid(grid, row, col, num):
    # Check row and column uniqueness
    for i in range(7):
        if grid[row][i] == num or grid[i][col] == num:
            return False
    
    # Check horizontal constraints
    if col > 0 and horizontal_constraints[row][col - 1] == '<' and grid[row][col - 1] != 0 and grid[row][col - 1] >= num:
        return False
    if col > 0 and horizontal_constraints[row][col - 1] == '>' and grid[row][col - 1] != 0 and grid[row][col - 1] <= num:
        return False
    if col < 6 and horizontal_constraints[row][col] == '<' and grid[row][col + 1] != 0 and grid[row][col + 1] <= num:
        return False
    if col < 6 and horizontal_constraints[row][col] == '>' and grid[row][col + 1] != 0 and grid[row][col + 1] >= num:
        return False
    
    # Check vertical constraints
    if row > 0 and vertical_constraints[row - 1][col] == '∧' and grid[row - 1][col] != 0 and grid[row - 1][col] >= num:
        return False
    if row > 0 and vertical_constraints[row - 1][col] == '∨' and grid[row - 1][col] != 0 and grid[row - 1][col] <= num:
        return False
    if row < 6 and vertical_constraints[row][col] == '∧' and grid[row + 1][col] != 0 and grid[row + 1][col] <= num:
        return False
    if row < 6 and vertical_constraints[row][col] == '∨' and grid[row + 1][col] != 0 and grid[row + 1][col] >= num:
        return False
    
    return True

def solve(grid):
    for row in range(7):
        for col in range(7):
            if grid[row][col] == 0:
                for num in range(1, 8):
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