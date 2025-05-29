def is_valid(grid, row, col, num):
    # Check if num is not in the current row and column
    for i in range(4):
        if grid[row][i] == num or grid[i][col] == num:
            return False
    
    # Check inequality constraints
    if row == 2 and col == 1 and grid[3][1] != 0 and not (num > grid[3][1]):
        return False
    if row == 3 and col == 1 and grid[2][1] != 0 and not (grid[2][1] > num):
        return False
    if row == 2 and col == 2 and grid[3][2] != 0 and not (num < grid[3][2]):
        return False
    if row == 3 and col == 2 and grid[2][2] != 0 and not (grid[2][2] < num):
        return False
    
    return True

def solve_futoshiki(grid):
    for row in range(4):
        for col in range(4):
            if grid[row][col] == 0:
                for num in range(1, 5):
                    if is_valid(grid, row, col, num):
                        grid[row][col] = num
                        if solve_futoshiki(grid):
                            return True
                        grid[row][col] = 0
                return False
    return True

# Initial grid with 0 representing empty cells
grid = [
    [4, 0, 0, 0],
    [3, 0, 0, 2],
    [0, 0, 0, 0],
    [0, 0, 3, 0]
]

solve_futoshiki(grid)
print(grid)