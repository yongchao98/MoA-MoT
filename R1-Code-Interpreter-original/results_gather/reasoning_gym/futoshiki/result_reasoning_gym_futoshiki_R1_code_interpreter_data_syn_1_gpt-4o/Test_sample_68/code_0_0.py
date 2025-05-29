from itertools import permutations

def is_valid(grid, row, col, num):
    # Check if num is not in the current row and column
    for i in range(5):
        if grid[row][i] == num or grid[i][col] == num:
            return False
    
    # Check horizontal inequalities
    if col > 0 and grid[row][col - 1] != 0:
        if grid[row][col - 1] < num and (row, col - 1, row, col) in inequalities:
            return False
        if grid[row][col - 1] > num and (row, col, row, col - 1) in inequalities:
            return False
    
    if col < 4 and grid[row][col + 1] != 0:
        if grid[row][col + 1] < num and (row, col + 1, row, col) in inequalities:
            return False
        if grid[row][col + 1] > num and (row, col, row, col + 1) in inequalities:
            return False
    
    # Check vertical inequalities
    if row > 0 and grid[row - 1][col] != 0:
        if grid[row - 1][col] < num and (row - 1, col, row, col) in inequalities:
            return False
        if grid[row - 1][col] > num and (row, col, row - 1, col) in inequalities:
            return False
    
    if row < 4 and grid[row + 1][col] != 0:
        if grid[row + 1][col] < num and (row + 1, col, row, col) in inequalities:
            return False
        if grid[row + 1][col] > num and (row, col, row + 1, col) in inequalities:
            return False
    
    return True

def solve_futoshiki(grid):
    for row in range(5):
        for col in range(5):
            if grid[row][col] == 0:
                for num in range(1, 6):
                    if is_valid(grid, row, col, num):
                        grid[row][col] = num
                        if solve_futoshiki(grid):
                            return True
                        grid[row][col] = 0
                return False
    return True

# Initial grid with 0 representing empty cells
grid = [
    [0, 0, 3, 2, 0],
    [0, 0, 0, 0, 0],
    [0, 0, 0, 1, 0],
    [0, 3, 0, 0, 0],
    [0, 0, 1, 0, 0]
]

# Inequalities (row1, col1, row2, col2) means grid[row1][col1] > grid[row2][col2]
inequalities = {
    (1, 2, 1, 3), (1, 3, 1, 4),
    (2, 0, 2, 1),
    (3, 0, 3, 1),
    (4, 3, 4, 4),
    (0, 0, 1, 0), (2, 0, 3, 0),
    (1, 4, 2, 4)
}

solve_futoshiki(grid)

for row in grid:
    print(row)