def is_valid(grid, row, col, num):
    # Check if num is not in the current row and column
    for i in range(4):
        if grid[row][i] == num or grid[i][col] == num:
            return False
    # Check inequality constraints
    if row == 1 and col == 0 and not (num > grid[2][0]):
        return False
    if row == 1 and col == 1 and not (num > grid[2][1]):
        return False
    return True

def solve_futoshiki(grid, row=0, col=0):
    if row == 4:
        return True
    if col == 4:
        return solve_futoshiki(grid, row + 1, 0)
    if grid[row][col] != 0:
        return solve_futoshiki(grid, row, col + 1)
    
    for num in range(1, 5):
        if is_valid(grid, row, col, num):
            grid[row][col] = num
            if solve_futoshiki(grid, row, col + 1):
                return True
            grid[row][col] = 0
    return False

# Initial grid with known values
grid = [
    [0, 0, 0, 0],
    [0, 0, 2, 0],
    [0, 4, 0, 0],
    [0, 0, 3, 2]
]

if solve_futoshiki(grid):
    for row in grid:
        print(row)
else:
    print("No solution found")