def is_valid(grid, row, col, num):
    # Check if num is not in the current row and column
    for i in range(9):
        if grid[row][i] == num or grid[i][col] == num:
            return False
    
    # Check horizontal constraints
    if col < 8:
        if grid[row][col + 1] != 0 and grid[row][col + 1] < num:
            return False
    if col > 0:
        if grid[row][col - 1] != 0 and grid[row][col - 1] > num:
            return False
    
    # Check vertical constraints
    if row < 8:
        if grid[row + 1][col] != 0 and grid[row + 1][col] < num:
            return False
    if row > 0:
        if grid[row - 1][col] != 0 and grid[row - 1][col] > num:
            return False
    
    return True

def solve_futoshiki(grid):
    for row in range(9):
        for col in range(9):
            if grid[row][col] == 0:
                for num in range(1, 10):
                    if is_valid(grid, row, col, num):
                        grid[row][col] = num
                        if solve_futoshiki(grid):
                            return True
                        grid[row][col] = 0
                return False
    return True

# Initial grid with constraints
grid = [
    [0, 9, 0, 0, 8, 0, 0, 3, 0],
    [0, 2, 0, 0, 9, 0, 0, 7, 0],
    [2, 0, 0, 8, 0, 0, 0, 4, 3],
    [4, 1, 0, 0, 0, 8, 0, 0, 0],
    [8, 4, 0, 5, 0, 0, 6, 2, 7],
    [0, 7, 0, 1, 5, 6, 0, 0, 4],
    [0, 0, 0, 0, 0, 1, 7, 9, 8],
    [6, 0, 0, 0, 0, 3, 4, 0, 0],
    [9, 6, 0, 7, 0, 2, 8, 5, 0]
]

solve_futoshiki(grid)

for row in grid:
    print(row)