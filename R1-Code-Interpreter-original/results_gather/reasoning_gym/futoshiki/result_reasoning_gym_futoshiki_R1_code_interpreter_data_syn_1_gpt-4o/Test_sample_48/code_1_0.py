def is_valid(grid, row, col, num):
    # Check row and column constraints
    for i in range(6):
        if grid[row][i] == num or grid[i][col] == num:
            return False
    
    # Check horizontal constraints
    if col > 0:
        if (row == 0 and col == 1 and grid[row][col-1] != 0 and grid[row][col-1] <= num) or \
           (row == 0 and col == 4 and grid[row][col-1] != 0 and grid[row][col-1] >= num) or \
           (row == 2 and col == 1 and grid[row][col-1] != 0 and grid[row][col-1] >= num) or \
           (row == 3 and col == 2 and grid[row][col-1] != 0 and grid[row][col-1] <= num) or \
           (row == 3 and col == 4 and grid[row][col-1] != 0 and grid[row][col-1] >= num) or \
           (row == 5 and col == 1 and grid[row][col-1] != 0 and grid[row][col-1] <= num):
            return False
    if col < 5:
        if (row == 0 and col == 0 and grid[row][col+1] != 0 and grid[row][col+1] >= num) or \
           (row == 0 and col == 3 and grid[row][col+1] != 0 and grid[row][col+1] <= num) or \
           (row == 2 and col == 0 and grid[row][col+1] != 0 and grid[row][col+1] <= num) or \
           (row == 2 and col == 3 and grid[row][col+1] != 0 and grid[row][col+1] >= num) or \
           (row == 3 and col == 1 and grid[row][col+1] != 0 and grid[row][col+1] >= num) or \
           (row == 3 and col == 3 and grid[row][col+1] != 0 and grid[row][col+1] <= num) or \
           (row == 5 and col == 0 and grid[row][col+1] != 0 and grid[row][col+1] >= num):
            return False
    
    # Check vertical constraints
    if row > 0:
        if (col == 0 and row == 1 and grid[row-1][col] != 0 and grid[row-1][col] <= num) or \
           (col == 2 and row == 1 and grid[row-1][col] != 0 and grid[row-1][col] <= num) or \
           (col == 5 and row == 1 and grid[row-1][col] != 0 and grid[row-1][col] <= num):
            return False
    if row < 5:
        if (col == 0 and row == 0 and grid[row+1][col] != 0 and grid[row+1][col] >= num) or \
           (col == 2 and row == 0 and grid[row+1][col] != 0 and grid[row+1][col] >= num) or \
           (col == 5 and row == 0 and grid[row+1][col] != 0 and grid[row+1][col] >= num):
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
    [0, 0, 0, 3, 0, 0],
    [0, 0, 6, 0, 1, 0],
    [0, 1, 0, 5, 0, 4],
    [0, 3, 0, 0, 6, 0],
    [4, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0]
]

solve_futoshiki(grid)

for row in grid:
    print(row)