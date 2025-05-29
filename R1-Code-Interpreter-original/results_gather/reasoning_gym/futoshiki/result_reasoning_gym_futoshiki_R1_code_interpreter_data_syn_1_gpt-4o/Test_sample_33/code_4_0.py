def is_valid(grid, row, col, num):
    # Check if num is not in the current row and column
    for i in range(7):
        if grid[row][i] == num or grid[i][col] == num:
            return False
    
    # Check horizontal constraints
    if col > 0:
        if (row == 0 and col == 1 and grid[row][col - 1] != 0 and grid[row][col - 1] >= num) or \
           (row == 0 and col == 2 and grid[row][col - 1] != 0 and grid[row][col - 1] <= num) or \
           (row == 2 and col == 5 and grid[row][col - 1] != 0 and grid[row][col - 1] >= num) or \
           (row == 4 and col == 1 and grid[row][col - 1] != 0 and grid[row][col - 1] <= num) or \
           (row == 5 and col == 4 and grid[row][col - 1] != 0 and grid[row][col - 1] >= num):
            return False
    
    # Check vertical constraints
    if row > 0:
        if (col == 3 and row == 1 and grid[row - 1][col] != 0 and grid[row - 1][col] >= num) or \
           (col == 6 and row == 2 and grid[row - 1][col] != 0 and grid[row - 1][col] >= num) or \
           (col == 6 and row == 5 and grid[row - 1][col] != 0 and grid[row - 1][col] <= num):
            return False
    
    return True

def solve_futoshiki(grid):
    for row in range(7):
        for col in range(7):
            if grid[row][col] == 0:
                for num in range(1, 8):
                    if is_valid(grid, row, col, num):
                        grid[row][col] = num
                        if solve_futoshiki(grid):
                            return True
                        grid[row][col] = 0
                return False
    return True

# Initial grid setup
grid = [
    [0, 0, 0, 0, 2, 0, 3],
    [0, 1, 0, 2, 0, 0, 0],
    [0, 0, 0, 6, 3, 0, 0],
    [0, 6, 4, 0, 0, 2, 0],
    [4, 0, 0, 0, 0, 0, 6],
    [0, 0, 0, 5, 0, 6, 1],
    [5, 0, 2, 0, 6, 0, 0]
]

if solve_futoshiki(grid):
    for row in grid:
        print(row)
else:
    print("No solution found")