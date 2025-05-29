def is_valid(grid, row, col, num):
    # Check if num is not in the current row
    if num in grid[row]:
        return False
    
    # Check if num is not in the current column
    for r in range(9):
        if grid[r][col] == num:
            return False
    
    # Check if num is not in the current 3x3 subgrid
    start_row, start_col = 3 * (row // 3), 3 * (col // 3)
    for r in range(start_row, start_row + 3):
        for c in range(start_col, start_col + 3):
            if grid[r][c] == num:
                return False
    
    return True

def solve_sudoku(grid):
    for row in range(9):
        for col in range(9):
            if grid[row][col] == 0:  # Find an empty cell
                for num in range(1, 10):  # Try numbers 1 to 9
                    if is_valid(grid, row, col, num):
                        grid[row][col] = num  # Place the number
                        if solve_sudoku(grid):  # Recursively solve the rest
                            return True
                        grid[row][col] = 0  # Backtrack
                return False  # No valid number found, trigger backtracking
    return True  # Solved

# Define the initial Sudoku grid
grid = [
    [0, 0, 0, 0, 4, 0, 6, 0, 9],
    [0, 6, 0, 5, 0, 0, 0, 0, 0],
    [0, 3, 0, 0, 0, 8, 0, 1, 0],
    [0, 0, 0, 0, 3, 0, 0, 0, 1],
    [0, 0, 0, 1, 0, 0, 0, 6, 7],
    [0, 0, 9, 0, 0, 2, 8, 0, 0],
    [0, 0, 0, 0, 0, 5, 0, 0, 2],
    [7, 0, 2, 0, 0, 0, 0, 0, 4],
    [8, 4, 0, 7, 0, 0, 0, 5, 6]
]

# Solve the Sudoku
solve_sudoku(grid)

# Print the solved Sudoku grid
for row in grid:
    print(row)