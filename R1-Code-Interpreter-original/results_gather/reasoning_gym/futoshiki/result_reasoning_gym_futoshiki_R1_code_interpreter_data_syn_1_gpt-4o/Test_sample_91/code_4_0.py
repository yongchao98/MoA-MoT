def is_valid_partial(grid, row, col, num):
    # Check if num is not in the current row
    if num in grid[row]:
        return False
    
    # Check if num is not in the current column
    for r in range(4):
        if grid[r][col] == num:
            return False
    
    # Check the constraints
    if row == 0 and col == 0 and num != 4:
        return False
    if row == 0 and col == 1 and grid[0][0] <= num:
        return False
    if row == 1 and col == 0 and grid[1][1] >= num:
        return False
    if row == 1 and col == 1 and (grid[1][0] <= num or grid[2][1] >= num):
        return False
    if row == 1 and col == 2 and grid[1][1] <= num:
        return False
    if row == 2 and col == 1 and grid[1][1] <= num:
        return False
    if row == 2 and col == 2 and num != 3:
        return False
    if row == 2 and col == 3 and grid[2][2] <= num:
        return False
    
    return True

def solve_futoshiki_backtrack(grid, row=0, col=0):
    if row == 4:
        return True
    
    next_row, next_col = (row, col + 1) if col < 3 else (row + 1, 0)
    
    if grid[row][col] != 0:
        return solve_futoshiki_backtrack(grid, next_row, next_col)
    
    for num in range(1, 5):
        if is_valid_partial(grid, row, col, num):
            grid[row][col] = num
            if solve_futoshiki_backtrack(grid, next_row, next_col):
                return True
            grid[row][col] = 0
    
    return False

# Initialize the grid with zeros
grid = [[0, 0, 0, 0] for _ in range(4)]

# Solve the puzzle
if solve_futoshiki_backtrack(grid):
    for row in grid:
        print(row)
else:
    print("No solution found")