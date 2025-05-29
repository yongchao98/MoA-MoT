def is_valid(grid, row, col, num):
    # Check if num is not in the current row and column
    for i in range(7):
        if grid[row][i] == num or grid[i][col] == num:
            return False
    
    # Check horizontal constraints
    if col > 0 and grid[row][col - 1] != 0:
        if grid[row][col - 1] < num and (row, col - 1, '<') in constraints:
            return False
        if grid[row][col - 1] > num and (row, col - 1, '>') in constraints:
            return False
    if col < 6 and grid[row][col + 1] != 0:
        if grid[row][col + 1] < num and (row, col, '<') in constraints:
            return False
        if grid[row][col + 1] > num and (row, col, '>') in constraints:
            return False
    
    # Check vertical constraints
    if row > 0 and grid[row - 1][col] != 0:
        if grid[row - 1][col] < num and (row - 1, col, '∧') in constraints:
            return False
        if grid[row - 1][col] > num and (row - 1, col, '∨') in constraints:
            return False
    if row < 6 and grid[row + 1][col] != 0:
        if grid[row + 1][col] < num and (row, col, '∧') in constraints:
            return False
        if grid[row + 1][col] > num and (row, col, '∨') in constraints:
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

# Initial grid setup
grid = [
    [3, 0, 5, 2, 0, 0, 0],
    [2, 0, 0, 6, 0, 0, 4],
    [0, 2, 0, 3, 0, 6, 7],
    [0, 4, 0, 0, 0, 0, 0],
    [6, 0, 0, 0, 5, 0, 0],
    [0, 0, 0, 0, 0, 1, 0],
    [0, 0, 4, 0, 3, 2, 0]
]

# Constraints as tuples (row, col, type)
constraints = {
    (0, 5, '<'), (1, 5, '<'), (3, 0, '<'), (3, 1, '<'), (3, 3, '<'), (4, 4, '<'),
    (5, 2, '<'), (6, 2, '>'), (1, 4, '∧'), (1, 6, '∧'), (2, 1, '∧'), (2, 2, '∧'),
    (4, 1, '∧'), (5, 1, '∧')
}

solve_futoshiki(grid)

# Print the solved grid
for row in grid:
    print(' '.join(str(num) if num != 0 else '_' for num in row))