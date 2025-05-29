def is_valid(grid, row, col, num):
    # Check if num is not in the current row
    if num in grid[row]:
        return False
    # Check if num is not in the current column
    for r in range(4):
        if grid[r][col] == num:
            return False
    return True

def solve_futoshiki(grid):
    for row in range(4):
        for col in range(4):
            if grid[row][col] == 0:  # Empty cell
                for num in range(1, 5):
                    if is_valid(grid, row, col, num):
                        grid[row][col] = num
                        if solve_futoshiki(grid):
                            return True
                        grid[row][col] = 0  # Backtrack
                return False
    return True

# Initial grid with 0 representing empty cells
grid = [
    [1, 0, 0, 0],
    [0, 0, 3, 0],
    [3, 0, 0, 0],
    [0, 1, 2, 0]
]

solve_futoshiki(grid)
print(grid)