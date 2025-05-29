def is_valid(grid, row, col, num):
    # Check if num is not in the current row
    if num in grid[row]:
        return False
    # Check if num is not in the current column
    for r in range(6):
        if grid[r][col] == num:
            return False
    return True

def solve_futoshiki(grid):
    for row in range(6):
        for col in range(6):
            if grid[row][col] == 0:  # 0 represents a blank cell
                for num in range(1, 7):
                    if is_valid(grid, row, col, num):
                        grid[row][col] = num
                        if solve_futoshiki(grid):
                            return True
                        grid[row][col] = 0  # Backtrack
                return False
    return True

# Initial grid with 0 representing blanks
grid = [
    [0, 2, 0, 0, 0, 4],
    [6, 1, 0, 0, 0, 5],
    [0, 6, 0, 2, 0, 0],
    [0, 0, 0, 0, 1, 0],
    [0, 0, 0, 3, 0, 2],
    [0, 4, 6, 0, 2, 0]
]

solve_futoshiki(grid)

# Print the solved grid
for row in grid:
    print(' '.join(map(str, row)))