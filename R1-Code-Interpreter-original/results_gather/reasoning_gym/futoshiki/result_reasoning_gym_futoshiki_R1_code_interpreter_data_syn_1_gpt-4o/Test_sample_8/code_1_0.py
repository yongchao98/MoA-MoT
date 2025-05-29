def is_valid(grid, row, col, num, constraints):
    # Check if num is not in the current row and column
    for i in range(7):
        if grid[row][i] == num or grid[i][col] == num:
            return False
    
    # Check horizontal constraints
    if col > 0 and grid[row][col - 1] != 0:
        if (row, col - 1, '<') in constraints and not (grid[row][col - 1] < num):
            return False
        if (row, col - 1, '>') in constraints and not (grid[row][col - 1] > num):
            return False
    if col < 6 and grid[row][col + 1] != 0:
        if (row, col, '<') in constraints and not (num < grid[row][col + 1]):
            return False
        if (row, col, '>') in constraints and not (num > grid[row][col + 1]):
            return False
    
    # Check vertical constraints
    if row > 0 and grid[row - 1][col] != 0:
        if (row - 1, col, '∧') in constraints and not (grid[row - 1][col] < num):
            return False
        if (row - 1, col, '∨') in constraints and not (grid[row - 1][col] > num):
            return False
    if row < 6 and grid[row + 1][col] != 0:
        if (row, col, '∧') in constraints and not (num < grid[row + 1][col]):
            return False
        if (row, col, '∨') in constraints and not (num > grid[row + 1][col]):
            return False
    
    return True

def solve_futoshiki(grid, constraints):
    for row in range(7):
        for col in range(7):
            if grid[row][col] == 0:
                for num in range(1, 8):
                    if is_valid(grid, row, col, num, constraints):
                        grid[row][col] = num
                        if solve_futoshiki(grid, constraints):
                            return True
                        grid[row][col] = 0
                return False
    return True

# Initial grid setup
grid = [
    [0, 6, 4, 0, 0, 0, 2],
    [7, 0, 0, 0, 1, 0, 0],
    [0, 2, 0, 0, 4, 6, 0],
    [4, 0, 2, 0, 0, 1, 0],
    [2, 0, 6, 0, 3, 0, 0],
    [0, 1, 3, 7, 0, 0, 0],
    [0, 0, 0, 3, 6, 0, 0]
]

# Constraints
constraints = {
    (0, 3, '<'),  # Row 1, Col 4 < Row 1, Col 5
    (4, 2, '∧')   # Row 5, Col 3 ∧ Row 6, Col 3
}

solve_futoshiki(grid, constraints)

# Print the solved grid
for row in grid:
    print(row)