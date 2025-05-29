from itertools import permutations

def is_valid(grid, row, col, num):
    # Check if num is not in the current row and column
    for i in range(6):
        if grid[row][i] == num or grid[i][col] == num:
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
    [0, 0, 6, 0, 1, 2],
    [6, 0, 0, 0, 3, 0],
    [0, 6, 0, 0, 0, 0],
    [0, 0, 0, 0, 2, 5],
    [1, 2, 3, 0, 0, 0],
    [0, 0, 2, 3, 5, 0]
]

# Solve the puzzle
solve_futoshiki(grid)

# Print the solved grid
for row in grid:
    print(row)