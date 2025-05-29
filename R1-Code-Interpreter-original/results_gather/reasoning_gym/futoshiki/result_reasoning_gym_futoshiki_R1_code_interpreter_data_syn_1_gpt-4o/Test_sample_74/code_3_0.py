def is_valid(grid, row, col, num, constraints):
    # Check if num is not in the current row and column
    for i in range(9):
        if grid[row][i] == num or grid[i][col] == num:
            return False
    
    # Check horizontal constraints
    if col > 0:
        if constraints[row][col - 1] == '<' and grid[row][col - 1] != 0 and grid[row][col - 1] >= num:
            return False
        if constraints[row][col - 1] == '>' and grid[row][col - 1] != 0 and grid[row][col - 1] <= num:
            return False
    if col < 8:
        if constraints[row][col] == '<' and grid[row][col + 1] != 0 and grid[row][col + 1] <= num:
            return False
        if constraints[row][col] == '>' and grid[row][col + 1] != 0 and grid[row][col + 1] >= num:
            return False
    
    # Check vertical constraints
    if row > 0:
        if constraints[row - 1][col] == '∧' and grid[row - 1][col] != 0 and grid[row - 1][col] >= num:
            return False
        if constraints[row - 1][col] == '∨' and grid[row - 1][col] != 0 and grid[row - 1][col] <= num:
            return False
    if row < 8:
        if constraints[row][col] == '∧' and grid[row + 1][col] != 0 and grid[row + 1][col] <= num:
            return False
        if constraints[row][col] == '∨' and grid[row + 1][col] != 0 and grid[row + 1][col] >= num:
            return False
    
    return True

def solve_futoshiki(grid, constraints):
    for row in range(9):
        for col in range(9):
            if grid[row][col] == 0:
                for num in range(1, 10):
                    if is_valid(grid, row, col, num, constraints):
                        grid[row][col] = num
                        if solve_futoshiki(grid, constraints):
                            return True
                        grid[row][col] = 0
                return False
    return True

# Initial grid setup with 0 representing empty cells
grid = [
    [0, 0, 4, 2, 0, 0, 0, 8, 0],
    [0, 4, 5, 3, 0, 0, 6, 0, 0],
    [0, 0, 9, 0, 0, 0, 0, 0, 0],
    [8, 0, 0, 0, 9, 7, 0, 0, 6],
    [0, 0, 0, 5, 6, 0, 0, 1, 2],
    [0, 0, 0, 0, 8, 6, 5, 0, 0],
    [0, 8, 0, 0, 0, 0, 4, 0, 7],
    [0, 0, 8, 0, 0, 0, 0, 1, 3],
    [6, 9, 0, 8, 0, 1, 0, 0, 4]
]

# Constraints setup
constraints = [
    ['<', '', '', '<', '', '<', '', '', ''],
    ['', '', '>', '', '', '', '', '<', ''],
    ['', '<', '>', '', '', '>', '<', '', ''],
    ['', '', '', '<', '', '', '', '', ''],
    ['', '', '', '', '', '', '>', '', ''],
    ['', '', '<', '', '', '', '>', '<', ''],
    ['', '>', '<', '>', '', '', '', '<', ''],
    ['', '<', '', '>', '', '', '>', '<', ''],
    ['<', '', '', '', '', '<', '', '>', '']
]

solve_futoshiki(grid, constraints)

# Print the solved grid
for row in grid:
    print(row)