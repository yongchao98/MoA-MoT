def is_valid(grid, row, col, num):
    # Check if num is not in the current row and column
    for i in range(8):
        if grid[row][i] == num or grid[i][col] == num:
            return False
    
    # Check horizontal constraints
    if col > 0 and grid[row][col - 1] != 0:
        if grid[row][col - 1] > num and (row, col - 1, '>') in constraints:
            return False
        if grid[row][col - 1] < num and (row, col - 1, '<') in constraints:
            return False
    if col < 7 and grid[row][col + 1] != 0:
        if grid[row][col + 1] < num and (row, col, '>') in constraints:
            return False
        if grid[row][col + 1] > num and (row, col, '<') in constraints:
            return False
    
    # Check vertical constraints
    if row > 0 and grid[row - 1][col] != 0:
        if grid[row - 1][col] > num and (row - 1, col, '∧') in constraints:
            return False
        if grid[row - 1][col] < num and (row - 1, col, '∨') in constraints:
            return False
    if row < 7 and grid[row + 1][col] != 0:
        if grid[row + 1][col] < num and (row, col, '∧') in constraints:
            return False
        if grid[row + 1][col] > num and (row, col, '∨') in constraints:
            return False
    
    return True

def solve_futoshiki(grid):
    for row in range(8):
        for col in range(8):
            if grid[row][col] == 0:
                for num in range(1, 9):
                    if is_valid(grid, row, col, num):
                        grid[row][col] = num
                        if solve_futoshiki(grid):
                            return True
                        grid[row][col] = 0
                return False
    return True

# Initial grid with 0 representing empty cells
grid = [
    [0, 7, 0, 2, 0, 0, 1, 5],
    [7, 0, 0, 0, 0, 4, 0, 2],
    [0, 5, 7, 0, 1, 0, 4, 6],
    [6, 0, 0, 3, 4, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 6, 0],
    [3, 2, 6, 0, 7, 0, 5, 0],
    [0, 4, 2, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 5, 0, 8, 1]
]

# Constraints in the form of (row, col, type)
constraints = {
    (0, 2, '>'), (1, 5, '∧'), (2, 3, '∨'), (3, 3, '<'), (4, 0, '<'),
    (4, 2, '∨'), (5, 3, '>'), (5, 4, '>'), (6, 0, '>'), (7, 4, '>')
}

solve_futoshiki(grid)

# Print the solved grid
for row in grid:
    print(row)