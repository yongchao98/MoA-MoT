def solve_sudoku(grid):
    """
    Solves a Sudoku puzzle using backtracking.
    """
    find = find_empty(grid)
    if not find:
        return True  # Puzzle is solved
    else:
        row, col = find

    for num in range(1, 10):
        if is_valid(grid, row, col, num):
            grid[row][col] = num

            if solve_sudoku(grid):
                return True

            # Backtrack if the solution is not found
            grid[row][col] = 0

    return False

def find_empty(grid):
    """
    Finds an empty cell (represented by 0) in the grid.
    """
    for i in range(9):
        for j in range(9):
            if grid[i][j] == 0:
                return (i, j)  # row, col
    return None

def is_valid(grid, row, col, num):
    """
    Checks if placing a number in a given cell is valid.
    """
    # Check row
    for j in range(9):
        if grid[row][j] == num and col != j:
            return False

    # Check column
    for i in range(9):
        if grid[i][col] == num and row != i:
            return False

    # Check 3x3 box
    box_x = col // 3
    box_y = row // 3
    for i in range(box_y * 3, box_y * 3 + 3):
        for j in range(box_x * 3, box_x * 3 + 3):
            if grid[i][j] == num and (i, j) != (row, col):
                return False

    return True

# Parse the bmatrix representation into a 9x9 grid
# 0 represents an empty cell
grid = [
    [5, 0, 0, 0, 8, 0, 0, 4, 9],
    [0, 0, 0, 5, 0, 0, 0, 3, 0],
    [0, 6, 7, 3, 0, 0, 0, 0, 1],
    [1, 5, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 2, 0, 8, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 1, 8],
    [7, 0, 0, 0, 0, 4, 1, 5, 0],
    [0, 3, 0, 0, 0, 2, 0, 0, 0],
    [4, 9, 0, 0, 5, 0, 0, 0, 3]
]

# Solve the puzzle
if solve_sudoku(grid):
    # Extract the top horizontal line from the solved grid
    top_line = grid[0]
    # Print the numbers separated by spaces
    print(' '.join(map(str, top_line)))
else:
    print("No solution exists for the given Sudoku.")
