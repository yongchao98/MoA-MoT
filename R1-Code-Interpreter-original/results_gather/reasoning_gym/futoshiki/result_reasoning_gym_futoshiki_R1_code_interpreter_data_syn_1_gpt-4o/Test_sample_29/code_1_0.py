def is_valid(grid, row, col, num):
    # Check if num is not in the current row and column
    for i in range(7):
        if grid[row][i] == num or grid[i][col] == num:
            return False
    
    # Check horizontal constraints
    if col > 0:
        if constraints[row][col - 1] == '<' and grid[row][col - 1] != 0 and grid[row][col - 1] >= num:
            return False
        if constraints[row][col - 1] == '>' and grid[row][col - 1] != 0 and grid[row][col - 1] <= num:
            return False
    if col < 6:
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
    if row < 6:
        if constraints[row][col] == '∧' and grid[row + 1][col] != 0 and grid[row + 1][col] <= num:
            return False
        if constraints[row][col] == '∨' and grid[row + 1][col] != 0 and grid[row + 1][col] >= num:
            return False
    
    return True

def solve_futoshiki(grid, row=0, col=0):
    if row == 7:
        return True
    if col == 7:
        return solve_futoshiki(grid, row + 1, 0)
    if grid[row][col] != 0:
        return solve_futoshiki(grid, row, col + 1)
    
    for num in range(1, 8):
        if is_valid(grid, row, col, num):
            grid[row][col] = num
            if solve_futoshiki(grid, row, col + 1):
                return True
            grid[row][col] = 0
    
    return False

# Initial grid with 0 representing empty cells
grid = [
    [0, 0, 0, 0, 3, 0, 0],
    [7, 5, 0, 0, 0, 0, 0],
    [0, 0, 0, 2, 0, 0, 3],
    [5, 7, 0, 0, 0, 1, 0],
    [0, 0, 7, 0, 0, 6, 2],
    [0, 0, 3, 0, 0, 0, 1],
    [3, 2, 0, 0, 5, 0, 0]
]

# Constraints for each cell
constraints = [
    ['', '∨', '∨', '', '', '', ''],
    ['', '', '', '<', '', '', '∨'],
    ['', '', '', '', '', '>', ''],
    ['', '>', '', '', '', '', ''],
    ['', '', '', '', '', '', ''],
    ['', '', '', '', '', '', ''],
    ['', '', '', '>', '', '', '']
]

solve_futoshiki(grid)

# Print the solved grid
for row in grid:
    print(row)