def is_valid(grid, row, col, num):
    # Check if num is not in the current row and column
    for i in range(8):
        if grid[row][i] == num or grid[i][col] == num:
            return False

    # Check horizontal constraints
    if col > 0 and grid[row][col - 1] != 0:
        if constraints[row][col - 1] == '<' and not (grid[row][col - 1] < num):
            return False
        if constraints[row][col - 1] == '>' and not (grid[row][col - 1] > num):
            return False

    if col < 7 and grid[row][col + 1] != 0:
        if constraints[row][col] == '<' and not (num < grid[row][col + 1]):
            return False
        if constraints[row][col] == '>' and not (num > grid[row][col + 1]):
            return False

    # Check vertical constraints
    if row > 0 and grid[row - 1][col] != 0:
        if constraints[row - 1][col] == '∧' and not (grid[row - 1][col] < num):
            return False
        if constraints[row - 1][col] == '∨' and not (grid[row - 1][col] > num):
            return False

    if row < 7 and grid[row + 1][col] != 0:
        if constraints[row][col] == '∧' and not (num < grid[row + 1][col]):
            return False
        if constraints[row][col] == '∨' and not (num > grid[row + 1][col]):
            return False

    return True

def solve_futoshiki(grid, row=0, col=0):
    if row == 8:
        return True
    if col == 8:
        return solve_futoshiki(grid, row + 1, 0)
    if grid[row][col] != 0:
        return solve_futoshiki(grid, row, col + 1)

    for num in range(1, 9):
        if is_valid(grid, row, col, num):
            grid[row][col] = num
            if solve_futoshiki(grid, row, col + 1):
                return True
            grid[row][col] = 0

    return False

# Initial grid with 0 representing empty cells
grid = [
    [0, 0, 0, 4, 0, 0, 0, 2],
    [2, 1, 8, 0, 0, 4, 3, 0],
    [8, 0, 2, 0, 0, 0, 1, 0],
    [3, 0, 0, 0, 1, 0, 6, 0],
    [6, 7, 0, 1, 0, 2, 0, 0],
    [0, 6, 3, 0, 7, 8, 0, 5],
    [7, 0, 0, 0, 3, 0, 0, 0],
    [4, 0, 0, 6, 0, 5, 0, 0]
]

# Constraints for the puzzle
constraints = [
    ['', '', '', '', '', '<', '', ''],
    ['∨', '', '', '', '', '', '', ''],
    ['', '', '', '', '', '', '', ''],
    ['', '', '', '', '', '', '>', ''],
    ['', '', '>', '', '', '', '', ''],
    ['', '<', '', '', '', '', '', ''],
    ['', '', '', '', '', '', '', ''],
    ['', '', '', '', '', '', '', '']
]

solve_futoshiki(grid)

for row in grid:
    print(row)