def is_valid(grid, row, col, num):
    # Check if num is not in the current row and column
    for i in range(6):
        if grid[row][i] == num or grid[i][col] == num:
            return False
    
    # Check horizontal inequalities
    if col > 0 and inequalities[row][col - 1] != '':
        if (inequalities[row][col - 1] == '>' and grid[row][col - 1] <= num) or \
           (inequalities[row][col - 1] == '<' and grid[row][col - 1] >= num):
            return False
    if col < 5 and inequalities[row][col] != '':
        if (inequalities[row][col] == '>' and num <= grid[row][col + 1]) or \
           (inequalities[row][col] == '<' and num >= grid[row][col + 1]):
            return False
    
    # Check vertical inequalities
    if row > 0 and vertical_inequalities[row - 1][col] != '':
        if (vertical_inequalities[row - 1][col] == '∧' and grid[row - 1][col] <= num) or \
           (vertical_inequalities[row - 1][col] == '∨' and grid[row - 1][col] >= num):
            return False
    if row < 5 and vertical_inequalities[row][col] != '':
        if (vertical_inequalities[row][col] == '∧' and num <= grid[row + 1][col]) or \
           (vertical_inequalities[row][col] == '∨' and num >= grid[row + 1][col]):
            return False
    
    return True

def solve_futoshiki(grid, row=0, col=0):
    if row == 6:
        return True
    if col == 6:
        return solve_futoshiki(grid, row + 1, 0)
    if grid[row][col] != 0:
        return solve_futoshiki(grid, row, col + 1)
    
    for num in range(1, 7):
        if is_valid(grid, row, col, num):
            grid[row][col] = num
            if solve_futoshiki(grid, row, col + 1):
                return True
            grid[row][col] = 0
    
    return False

# Initial grid with 0 representing empty cells
grid = [
    [0, 0, 0, 0, 0, 0],
    [0, 0, 4, 0, 0, 3],
    [4, 0, 5, 0, 0, 0],
    [0, 0, 3, 0, 4, 0],
    [0, 2, 0, 1, 0, 0],
    [0, 0, 0, 0, 0, 1]
]

# Horizontal inequalities
inequalities = [
    ['>', '', '', '', '', ''],
    ['<', '', '>', '', '', ''],
    ['', '>', '', '>', '', ''],
    ['', '<', '', '', '<', ''],
    ['', '<', '', '', '', ''],
    ['', '', '<', '', '', '']
]

# Vertical inequalities
vertical_inequalities = [
    ['', '', '', '', '', ''],
    ['∧', '', '', '', '∨', '∧'],
    ['', '', '', '', '', '∧'],
    ['', '', '∧', '', '', ''],
    ['', '∧', '', '∨', '', ''],
    ['', '', '', '', '', '']
]

solve_futoshiki(grid)

# Print the solved grid
for row in grid:
    print(' '.join(map(str, row)))