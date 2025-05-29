def is_valid(grid, row, col, num):
    # Check if num is not in the current row and column
    for i in range(9):
        if grid[row][i] == num or grid[i][col] == num:
            return False

    # Check horizontal constraints
    if col > 0 and grid[row][col - 1] != 0:
        if (grid[row][col - 1] < num and '<' in constraints[row][col - 1]) or \
           (grid[row][col - 1] > num and '>' in constraints[row][col - 1]):
            return False
    if col < 8 and grid[row][col + 1] != 0:
        if (grid[row][col + 1] < num and '>' in constraints[row][col]) or \
           (grid[row][col + 1] > num and '<' in constraints[row][col]):
            return False

    # Check vertical constraints
    if row > 0 and grid[row - 1][col] != 0:
        if (grid[row - 1][col] < num and '∧' in constraints[row - 1][col]) or \
           (grid[row - 1][col] > num and '∨' in constraints[row - 1][col]):
            return False
    if row < 8 and grid[row + 1][col] != 0:
        if (grid[row + 1][col] < num and '∨' in constraints[row][col]) or \
           (grid[row + 1][col] > num and '∧' in constraints[row][col]):
            return False

    return True

def solve_futoshiki(grid, row=0, col=0):
    if row == 9:
        return True
    if col == 9:
        return solve_futoshiki(grid, row + 1, 0)
    if grid[row][col] != 0:
        return solve_futoshiki(grid, row, col + 1)

    for num in range(1, 10):
        if is_valid(grid, row, col, num):
            grid[row][col] = num
            if solve_futoshiki(grid, row, col + 1):
                return True
            grid[row][col] = 0

    return False

# Initial grid setup with 0 for blanks
grid = [
    [7, 0, 0, 8, 9, 0, 0, 6, 5, 0, 0],
    [5, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 9, 7, 2, 0, 3],
    [3, 0, 1, 0, 0, 4, 0, 0, 0],
    [2, 1, 0, 3, 0, 0, 0, 9, 0, 6],
    [9, 7, 6, 0, 0, 0, 0, 3, 0],
    [6, 0, 5, 4, 7, 0, 3, 0, 1, 0],
    [0, 5, 0, 0, 4, 0, 0, 7, 0],
    [0, 0, 0, 5, 0, 0, 0, 0, 0, 8]
]

# Constraints for the puzzle
constraints = [
    ['', '<', '', '>', '', '<', '', '>', '', '<', ''],
    ['>', '', '<', '', '<', '', '', '', ''],
    ['<', '', '', '', '', '', '', '', ''],
    ['', '', '<', '', '', '', '<', '', ''],
    ['', '<', '', '', '', '<', '', '', ''],
    ['>', '', '', '', '', '', '>', '', ''],
    ['', '', '', '>', '', '<', '', '>', '', ''],
    ['', '>', '', '', '', '', '', '', ''],
    ['', '', '', '', '', '', '>', '', '<']
]

solve_futoshiki(grid)

# Print the solved grid
for row in grid:
    print(row)