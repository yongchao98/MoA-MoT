def is_valid(grid, row, col, num):
    # Check if num is not in the current row
    if num in grid[row]:
        return False
    
    # Check if num is not in the current column
    for r in range(5):
        if grid[r][col] == num:
            return False
    
    # Check the constraints
    if row == 0 and col == 1 and grid[1][1] != 0 and not (num < grid[1][1]):
        return False
    if row == 1 and col == 1 and grid[0][1] != 0 and not (grid[0][1] < num):
        return False
    if row == 1 and col == 0 and grid[2][0] != 0 and not (num > grid[2][0]):
        return False
    if row == 2 and col == 0 and grid[1][0] != 0 and not (grid[1][0] > num):
        return False
    if row == 2 and col == 0 and num <= 2:
        return False
    
    return True

def solve_futoshiki(grid, row=0, col=0):
    if row == 5:
        return True
    
    next_row, next_col = (row, col + 1) if col < 4 else (row + 1, 0)
    
    if grid[row][col] != 0:
        return solve_futoshiki(grid, next_row, next_col)
    
    for num in range(1, 6):
        if is_valid(grid, row, col, num):
            grid[row][col] = num
            if solve_futoshiki(grid, next_row, next_col):
                return True
            grid[row][col] = 0
    
    return False

# Initial grid with known values
grid = [
    [0, 3, 0, 5, 0],
    [0, 0, 0, 0, 0],
    [0, 0, 5, 1, 0],
    [0, 4, 2, 0, 0],
    [3, 0, 0, 0, 0]
]

solve_futoshiki(grid)
for row in grid:
    print(row)