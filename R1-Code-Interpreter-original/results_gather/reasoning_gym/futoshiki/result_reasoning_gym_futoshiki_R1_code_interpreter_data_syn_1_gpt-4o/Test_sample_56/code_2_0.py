def is_valid(grid, row, col, num):
    # Check if num is not in the current row and column
    for i in range(7):
        if grid[row][i] == num or grid[i][col] == num:
            return False
    
    # Check horizontal constraints
    if col > 0 and grid[row][col - 1] != 0:
        if (row, col - 1, '>') in constraints and grid[row][col - 1] <= num:
            return False
        if (row, col - 1, '<') in constraints and grid[row][col - 1] >= num:
            return False
    if col < 6 and grid[row][col + 1] != 0:
        if (row, col, '<') in constraints and grid[row][col + 1] <= num:
            return False
        if (row, col, '>') in constraints and grid[row][col + 1] >= num:
            return False
    
    # Check vertical constraints
    if row > 0 and grid[row - 1][col] != 0:
        if (row - 1, col, '∨') in constraints and grid[row - 1][col] <= num:
            return False
        if (row - 1, col, '∧') in constraints and grid[row - 1][col] >= num:
            return False
    if row < 6 and grid[row + 1][col] != 0:
        if (row, col, '∧') in constraints and grid[row + 1][col] <= num:
            return False
        if (row, col, '∨') in constraints and grid[row + 1][col] >= num:
            return False
    
    return True

def solve_futoshiki(grid):
    for row in range(7):
        for col in range(7):
            if grid[row][col] == 0:
                for num in range(1, 8):
                    if is_valid(grid, row, col, num):
                        grid[row][col] = num
                        if solve_futoshiki(grid):
                            return True
                        grid[row][col] = 0
                return False
    return True

# Initial grid with 0 representing empty cells
grid = [
    [4, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 6, 0, 1, 0],
    [0, 0, 0, 0, 0, 7, 4],
    [0, 0, 5, 0, 2, 3, 0],
    [6, 0, 0, 1, 0, 0, 7],
    [0, 2, 0, 7, 0, 5, 0],
    [3, 0, 7, 0, 4, 0, 0]
]

# Constraints in the form (row, col, type)
constraints = [
    (0, 3, '>'), (6, 0, '>'),
    (1, 3, '∧'), (3, 4, '∧'),
    (1, 5, '∨'), (2, 6, '∨'), (3, 6, '∨')
]

solve_futoshiki(grid)

# Print the solved grid
for row in grid:
    print(' '.join(map(str, row)))