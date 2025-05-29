# Define the initial Sudoku grid
sudoku_grid = [
    [0, 0, 6, 0, 0, 0, 5, 0, 4],
    [0, 0, 0, 3, 0, 0, 8, 0, 7],
    [0, 3, 4, 0, 0, 9, 0, 1, 0],
    [1, 8, 0, 2, 9, 0, 0, 0, 0],
    [3, 0, 0, 5, 0, 0, 9, 0, 0],
    [0, 6, 0, 0, 0, 0, 0, 0, 2],
    [7, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 3, 0, 0, 6, 9],
    [0, 0, 0, 9, 0, 0, 0, 8, 5]
]

# Function to check if a number can be placed in a given position
def is_valid(grid, row, col, num):
    # Check row
    for i in range(9):
        if grid[row][i] == num:
            return False
    # Check column
    for i in range(9):
        if grid[i][col] == num:
            return False
    # Check 3x3 subgrid
    start_row, start_col = 3 * (row // 3), 3 * (col // 3)
    for i in range(3):
        for j in range(3):
            if grid[start_row + i][start_col + j] == num:
                return False
    return True

# Function to solve the Sudoku using backtracking
def solve_sudoku(grid):
    for row in range(9):
        for col in range(9):
            if grid[row][col] == 0:  # Find an empty cell
                for num in range(1, 10):  # Try numbers 1 to 9
                    if is_valid(grid, row, col, num):
                        grid[row][col] = num  # Place the number
                        if solve_sudoku(grid):  # Recursively solve the rest
                            return True
                        grid[row][col] = 0  # Backtrack if needed
                return False  # Trigger backtracking
    return True

# Solve the Sudoku puzzle
solve_sudoku(sudoku_grid)

# Print the solved Sudoku grid
print(sudoku_grid)