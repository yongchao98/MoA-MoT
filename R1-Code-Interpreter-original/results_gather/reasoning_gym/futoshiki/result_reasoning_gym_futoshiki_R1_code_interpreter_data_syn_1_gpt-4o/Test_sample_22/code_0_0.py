def is_valid(grid, row, col, num):
    # Check if num is not in the current row and column
    for i in range(5):
        if grid[row][i] == num or grid[i][col] == num:
            return False

    # Check horizontal constraints
    if col > 0 and grid[row][col - 1] != 0 and grid[row][col - 1] > num:
        return False
    if col < 4 and grid[row][col + 1] != 0 and grid[row][col + 1] < num:
        return False

    # Check vertical constraints
    if row > 0 and grid[row - 1][col] != 0 and grid[row - 1][col] > num:
        return False
    if row < 4 and grid[row + 1][col] != 0 and grid[row + 1][col] < num:
        return False

    return True

def solve_futoshiki(grid):
    for row in range(5):
        for col in range(5):
            if grid[row][col] == 0:
                for num in range(1, 6):
                    if is_valid(grid, row, col, num):
                        grid[row][col] = num
                        if solve_futoshiki(grid):
                            return True
                        grid[row][col] = 0
                return False
    return True

# Initial grid with constraints
grid = [
    [0, 4, 2, 0, 0],
    [3, 0, 0, 0, 0],
    [0, 0, 0, 0, 0],
    [0, 2, 0, 0, 0],
    [0, 5, 0, 0, 0]
]

# Solve the puzzle
solve_futoshiki(grid)

# Print the solved grid
for row in grid:
    print(row)