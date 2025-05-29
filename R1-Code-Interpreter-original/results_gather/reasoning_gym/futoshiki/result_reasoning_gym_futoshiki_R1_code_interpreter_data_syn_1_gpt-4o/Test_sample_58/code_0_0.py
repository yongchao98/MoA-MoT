from itertools import permutations

def is_valid(grid, row, col, num):
    # Check if num is not in the current row and column
    for i in range(5):
        if grid[row][i] == num or grid[i][col] == num:
            return False
    
    # Check inequality constraints
    if row == 1 and col == 0 and not (num > grid[2][0]):
        return False
    if row == 2 and col == 3 and not (num > grid[3][3]):
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

# Initial grid with given numbers and 0 for blanks
grid = [
    [0, 0, 0, 3, 4],
    [5, 0, 0, 4, 0],
    [0, 0, 3, 0, 0],
    [0, 5, 0, 0, 1],
    [0, 0, 0, 0, 2]
]

solve_futoshiki(grid)
print(grid)