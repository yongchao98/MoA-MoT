def is_valid(grid, row, col, num):
    # Check if num is not in the current row and column
    for i in range(6):
        if grid[row][i] == num or grid[i][col] == num:
            return False
    
    # Check inequality constraints
    if row == 0 and col == 2 and grid[row][col + 1] != 0 and not (num > grid[row][col + 1]):
        return False
    if row == 0 and col == 3 and grid[row][col - 1] != 0 and not (grid[row][col - 1] > num):
        return False
    if row == 2 and col == 2 and grid[row + 1][col] != 0 and not (num < grid[row + 1][col]):
        return False
    if row == 3 and col == 2 and grid[row - 1][col] != 0 and not (grid[row - 1][col] < num):
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
    [0, 0, 0, 0, 3, 1],
    [0, 0, 5, 0, 0, 0],
    [0, 0, 2, 0, 0, 0],
    [5, 0, 3, 2, 6, 0],
    [0, 0, 0, 5, 1, 0],
    [4, 0, 0, 0, 0, 0]
]

solve_futoshiki(grid)
for row in grid:
    print(row)