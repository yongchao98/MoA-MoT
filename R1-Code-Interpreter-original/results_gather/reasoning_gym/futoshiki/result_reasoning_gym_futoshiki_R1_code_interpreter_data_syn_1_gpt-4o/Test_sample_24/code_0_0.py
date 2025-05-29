def is_valid(grid, row, col, num):
    # Check if num is not in the current row and column
    for i in range(8):
        if grid[row][i] == num or grid[i][col] == num:
            return False

    # Check the inequality constraints
    if row == 1 and col == 1 and grid[row][col + 1] != 0 and num <= grid[row][col + 1]:
        return False
    if row == 2 and col == 7 and grid[row + 1][col] != 0 and num <= grid[row + 1][col]:
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

# Initial grid setup
grid = [
    [0, 0, 0, 4, 0, 5, 6, 8],
    [6, 0, 0, 0, 0, 3, 5, 2],
    [8, 0, 3, 0, 7, 1, 4, 0],
    [0, 0, 8, 1, 4, 0, 0, 0],
    [0, 0, 1, 0, 5, 0, 0, 7],
    [4, 0, 0, 0, 0, 0, 8, 0],
    [0, 7, 2, 0, 0, 0, 0, 0],
    [0, 2, 0, 8, 0, 0, 0, 0]
]

solve_futoshiki(grid)

# Print the solved grid
for row in grid:
    print(row)