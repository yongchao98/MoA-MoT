def is_valid(grid, row, col, num):
    # Check if num is not in the current row and column
    for i in range(7):
        if grid[row][i] == num or grid[i][col] == num:
            return False
    
    # Check inequality constraints
    if (row, col) == (2, 5) and grid[3][5] != 0 and not (num > grid[3][5]):
        return False
    if (row, col) == (3, 5) and grid[2][5] != 0 and not (grid[2][5] > num):
        return False
    if (row, col) == (4, 2) and grid[5][2] != 0 and not (num > grid[5][2]):
        return False
    if (row, col) == (5, 2) and grid[4][2] != 0 and not (grid[4][2] > num):
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
    [0, 0, 0, 0, 7, 0, 0],
    [3, 2, 0, 0, 0, 0, 0],
    [0, 0, 2, 5, 0, 0, 7],
    [0, 0, 6, 0, 1, 5, 0],
    [4, 0, 0, 2, 6, 1, 0],
    [7, 0, 0, 0, 0, 3, 0],
    [0, 5, 0, 0, 0, 2, 4]
]

solve_futoshiki(grid)

for row in grid:
    print(row)