def is_valid(grid, row, col, num):
    # Check if num is not in the current row and column
    for i in range(6):
        if grid[row][i] == num or grid[i][col] == num:
            return False

    # Check the constraints
    if row == 3 and col == 4 and grid[row][5] != 0 and not (num > grid[row][5]):
        return False
    if row == 3 and col == 3 and grid[row][4] != 0 and not (num > grid[row][4]):
        return False

    return True

def solve_futoshiki(grid):
    for row in range(6):
        for col in range(6):
            if grid[row][col] == 0:
                for num in range(1, 7):
                    if is_valid(grid, row, col, num):
                        grid[row][col] = num
                        if solve_futoshiki(grid):
                            return True
                        grid[row][col] = 0
                return False
    return True

# Initial grid with 0 representing empty cells
grid = [
    [0, 0, 0, 0, 0, 3],
    [0, 0, 1, 0, 0, 0],
    [0, 0, 6, 3, 0, 0],
    [0, 2, 0, 0, 0, 0],
    [1, 0, 4, 0, 0, 5],
    [4, 0, 0, 0, 6, 0]
]

solve_futoshiki(grid)

# Print the solved grid
for row in grid:
    print(' '.join(map(str, row)))