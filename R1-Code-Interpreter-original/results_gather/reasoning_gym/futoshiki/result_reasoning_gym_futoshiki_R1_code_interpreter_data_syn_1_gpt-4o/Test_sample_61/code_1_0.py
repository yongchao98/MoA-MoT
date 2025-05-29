def is_valid(grid, row, col, num, row_constraints, col_constraints):
    # Check if num is already in the row or column
    for i in range(7):
        if grid[row][i] == num or grid[i][col] == num:
            return False
    
    # Check horizontal constraints
    if col > 0:
        if row_constraints[row][col - 1] == '>' and grid[row][col - 1] != 0 and grid[row][col - 1] <= num:
            return False
        if row_constraints[row][col - 1] == '<' and grid[row][col - 1] != 0 and grid[row][col - 1] >= num:
            return False
    if col < 6:
        if row_constraints[row][col] == '>' and grid[row][col + 1] != 0 and grid[row][col + 1] >= num:
            return False
        if row_constraints[row][col] == '<' and grid[row][col + 1] != 0 and grid[row][col + 1] <= num:
            return False
    
    # Check vertical constraints
    if row > 0:
        if col_constraints[row - 1][col] == '∧' and grid[row - 1][col] != 0 and grid[row - 1][col] <= num:
            return False
        if col_constraints[row - 1][col] == '∨' and grid[row - 1][col] != 0 and grid[row - 1][col] >= num:
            return False
    if row < 6:
        if col_constraints[row][col] == '∧' and grid[row + 1][col] != 0 and grid[row + 1][col] >= num:
            return False
        if col_constraints[row][col] == '∨' and grid[row + 1][col] != 0 and grid[row + 1][col] <= num:
            return False
    
    return True

def solve_futoshiki(grid, row_constraints, col_constraints):
    for row in range(7):
        for col in range(7):
            if grid[row][col] == 0:
                for num in range(1, 8):
                    if is_valid(grid, row, col, num, row_constraints, col_constraints):
                        grid[row][col] = num
                        if solve_futoshiki(grid, row_constraints, col_constraints):
                            return True
                        grid[row][col] = 0
                return False
    return True

# Initial grid with 0 representing empty cells
grid = [
    [0, 0, 0, 2, 0, 0, 0],
    [0, 0, 0, 0, 5, 0, 0],
    [0, 7, 0, 0, 0, 0, 0],
    [0, 0, 0, 6, 4, 2, 7],
    [4, 0, 0, 3, 0, 5, 0],
    [3, 0, 5, 0, 0, 0, 0],
    [0, 6, 0, 0, 0, 0, 3]
]

# Row constraints (horizontal)
row_constraints = [
    ['', '', '', '', '', '', ''],
    ['>', '', '', '>', '>', '<', ''],
    ['', '', '', '', '', '', ''],
    ['<', '', '', '', '', '', ''],
    ['<', '', '', '', '', '<', ''],
    ['', '<', '', '', '>', '', ''],
    ['', '', '', '', '', '', '']
]

# Column constraints (vertical)
col_constraints = [
    ['', '∨', '', '∧', '∧', '∨', ''],
    ['', '', '', '', '', '', ''],
    ['', '∨', '', '', '', '', ''],
    ['', '', '', '', '', '', ''],
    ['', '', '', '', '', '', '∨'],
    ['', '', '∨', '', '', '', ''],
    ['∧', '', '∨', '', '', '', '']
]

# Solve the puzzle
solve_futoshiki(grid, row_constraints, col_constraints)

# Print the solution
for row in grid:
    print(row)