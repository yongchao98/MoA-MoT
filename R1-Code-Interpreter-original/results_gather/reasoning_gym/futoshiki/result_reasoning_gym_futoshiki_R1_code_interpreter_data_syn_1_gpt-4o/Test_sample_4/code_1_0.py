def is_valid_partial(grid, row, col, num):
    # Check if num is not in the current row
    if num in grid[row]:
        return False
    # Check if num is not in the current column
    for r in range(6):
        if grid[r][col] == num:
            return False
    # Check the constraints
    if row == 0 and col == 1 and grid[0][2] != 0 and not (num > grid[0][2]):
        return False
    if row == 0 and col == 2 and grid[0][1] != 0 and not (grid[0][1] > num):
        return False
    if row == 2 and col == 2 and grid[3][2] != 0 and not (num > grid[3][2]):
        return False
    if row == 3 and col == 2 and grid[2][2] != 0 and not (grid[2][2] > num):
        return False
    return True

def solve_futoshiki_backtrack(grid, row=0, col=0):
    if row == 6:
        return True
    if col == 6:
        return solve_futoshiki_backtrack(grid, row + 1, 0)
    if grid[row][col] != 0:
        return solve_futoshiki_backtrack(grid, row, col + 1)
    
    for num in range(1, 7):
        if is_valid_partial(grid, row, col, num):
            grid[row][col] = num
            if solve_futoshiki_backtrack(grid, row, col + 1):
                return True
            grid[row][col] = 0
    return False

def solve_futoshiki():
    grid = [
        [0, 0, 0, 0, 3, 1],
        [0, 0, 5, 0, 0, 0],
        [0, 0, 2, 0, 0, 0],
        [5, 0, 3, 2, 6, 0],
        [0, 0, 0, 5, 1, 0],
        [4, 0, 0, 0, 0, 0]
    ]
    
    if solve_futoshiki_backtrack(grid):
        return grid
    else:
        return None

solution = solve_futoshiki()
for row in solution:
    print(row)