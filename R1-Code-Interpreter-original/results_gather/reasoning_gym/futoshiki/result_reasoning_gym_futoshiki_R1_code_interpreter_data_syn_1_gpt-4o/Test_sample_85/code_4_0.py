def is_valid(grid, row, col, num):
    # Check if num is not in the current row and column
    for i in range(5):
        if grid[row][i] == num or grid[i][col] == num:
            return False

    # Check specific constraints
    if row == 1 and col == 4 and num != 5:
        return False
    if row == 2 and col == 2 and num != 5:
        return False
    if row == 3:
        if col == 1 and num != 2:
            return False
        if col == 2 and num != 3:
            return False
        if col == 3 and num != 1:
            return False
    if row == 4 and col == 2 and num <= 2:
        return False

    # Check vertical constraints
    if row == 1 and col == 2 and grid[2][2] != 0 and num <= grid[2][2]:
        return False
    if row == 1 and col == 3 and grid[2][3] != 0 and num >= grid[2][3]:
        return False
    if row == 3 and col == 0 and grid[4][0] != 0 and num >= grid[4][0]:
        return False
    if row == 3 and col == 2 and grid[4][2] != 0 and num >= grid[4][2]:
        return False

    return True

def solve_futoshiki(grid, row=0, col=0):
    if row == 5:
        return True
    if col == 5:
        return solve_futoshiki(grid, row + 1, 0)
    if grid[row][col] != 0:
        return solve_futoshiki(grid, row, col + 1)

    for num in range(1, 6):
        if is_valid(grid, row, col, num):
            grid[row][col] = num
            if solve_futoshiki(grid, row, col + 1):
                return True
            grid[row][col] = 0

    return False

# Initial grid with 0 representing empty cells
grid = [
    [0, 0, 0, 0, 0],
    [0, 0, 0, 0, 5],
    [0, 0, 5, 0, 0],
    [0, 2, 3, 1, 0],
    [0, 0, 0, 2, 3]
]

if solve_futoshiki(grid):
    for row in grid:
        print(row)
else:
    print("No solution found")