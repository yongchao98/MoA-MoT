def is_valid(grid, row, col, num):
    # Check if num is not in the current row and column
    for i in range(5):
        if grid[row][i] == num or grid[i][col] == num:
            return False
    
    # Check the inequality constraint for the last row
    if row == 4 and col == 4 and not (num > grid[4][3]):
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

# Initial grid with known numbers
grid = [
    [0, 1, 0, 0, 3],
    [0, 0, 0, 3, 0],
    [0, 2, 0, 0, 5],
    [3, 0, 0, 1, 0],
    [0, 0, 4, 5, 0]
]

if solve_futoshiki(grid):
    for row in grid:
        print(row)
else:
    print("No solution found")