def is_valid(grid, row, col, num):
    # Check if num is not in the current row and column
    for i in range(7):
        if grid[row][i] == num or grid[i][col] == num:
            return False

    # Check horizontal inequalities
    if col > 0 and grid[row][col - 1] != 0:
        if (grid[row][col - 1] < num and inequalities[row][col - 1] == '>') or \
           (grid[row][col - 1] > num and inequalities[row][col - 1] == '<'):
            return False
    if col < 6 and grid[row][col + 1] != 0:
        if (grid[row][col + 1] < num and inequalities[row][col] == '<') or \
           (grid[row][col + 1] > num and inequalities[row][col] == '>'):
            return False

    # Check vertical inequalities
    if row > 0 and grid[row - 1][col] != 0:
        if (grid[row - 1][col] < num and inequalities[row - 1][col] == '∨') or \
           (grid[row - 1][col] > num and inequalities[row - 1][col] == '∧'):
            return False
    if row < 6 and grid[row + 1][col] != 0:
        if (grid[row + 1][col] < num and inequalities[row][col] == '∧') or \
           (grid[row + 1][col] > num and inequalities[row][col] == '∨'):
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

# Initial grid with given numbers and 0 for blanks
grid = [
    [4, 0, 0, 0, 0, 0, 2],
    [0, 6, 0, 0, 4, 0, 1, 0],
    [0, 0, 3, 0, 6, 0, 5],
    [0, 0, 0, 1, 0, 7, 3],
    [0, 0, 4, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 0],
    [0, 4, 0, 5, 2, 0, 0]
]

# Inequalities for the grid
inequalities = [
    ['', '∧', '∧', '', '', '∨', '∧'],
    ['', '>', '>', '>', '', '', ''],
    ['', '', '', '<', '', '', '∨'],
    ['>', '>', '>', '', '', '', ''],
    ['∨', '∧', '', '', '∧', '∨', '∧'],
    ['', '>', '', '', '∨', '', ''],
    ['', '<', '', '', '', '', '']
]

solve_futoshiki(grid)

# Print the solved grid
for row in grid:
    print(' '.join(map(str, row)))