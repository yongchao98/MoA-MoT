def solve_sudoku(grid):
    """
    Solves a Sudoku puzzle using backtracking.
    """
    find = find_empty(grid)
    if not find:
        return True  # Puzzle is solved if no empty cells are left
    else:
        row, col = find

    for num in range(1, 10):
        if is_valid(grid, num, (row, col)):
            grid[row][col] = num

            if solve_sudoku(grid):
                return True

            # Backtrack if the solution is not found
            grid[row][col] = 0

    return False

def find_empty(grid):
    """
    Finds an empty cell in the grid (represented by 0).
    """
    for i in range(9):
        for j in range(9):
            if grid[i][j] == 0:
                return (i, j)  # Returns (row, col)
    return None

def is_valid(grid, num, pos):
    """
    Checks if placing a number in a given position is valid.
    """
    row, col = pos

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
            if grid[i][j] == num and (i, j) != pos:
                return False

    return True

# Parse the input grid from the problem description
# 0 represents an empty cell
puzzle = [
    [5, 0, 0, 0, 8, 0, 0, 4, 9],
    [0, 0, 0, 5, 0, 0, 0, 3, 0],
    [0, 0, 6, 0, 7, 3, 0, 0, 1],
    [1, 0, 5, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 2, 0, 8, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 1, 8],
    [7, 0, 0, 0, 0, 4, 1, 5, 0],
    [0, 0, 3, 0, 0, 2, 0, 0, 0],
    [4, 0, 9, 0, 0, 5, 0, 0, 3]
]

# Solve the puzzle
if solve_sudoku(puzzle):
    # The puzzle is modified in-place.
    # Print the top horizontal line of the solved sudoku.
    # The format is 9 space-separated numbers.
    top_row = puzzle[0]
    # The final print statement includes each number from the solved row.
    print(f"{top_row[0]} {top_row[1]} {top_row[2]} {top_row[3]} {top_row[4]} {top_row[5]} {top_row[6]} {top_row[7]} {top_row[8]}")
else:
    print("No solution exists for the given Sudoku.")
