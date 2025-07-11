import sys

def solve_sudoku(grid):
    """
    Solves a Sudoku puzzle using a backtracking algorithm.
    The function modifies the grid in-place.
    """
    # Find the next empty cell (represented by 0)
    find = find_empty(grid)
    if not find:
        # If no empty cell is found, the puzzle is solved
        return True
    else:
        row, col = find

    # Try numbers from 1 to 9 in the empty cell
    for num in range(1, 10):
        # Check if the number is valid in the current position
        if is_valid(grid, num, (row, col)):
            # If valid, place the number
            grid[row][col] = num

            # Recursively call solve_sudoku for the next empty cell
            if solve_sudoku(grid):
                return True

            # If the recursive call fails, backtrack by resetting the cell
            grid[row][col] = 0

    # If no number from 1-9 works, return False to trigger backtracking
    return False

def find_empty(grid):
    """
    Finds the first empty cell (value 0) in the grid and returns its (row, col) position.
    Returns None if the grid is full.
    """
    for i in range(9):
        for j in range(9):
            if grid[i][j] == 0:
                return (i, j)
    return None

def is_valid(grid, num, pos):
    """
    Checks if placing 'num' at 'pos' (row, col) is valid according to Sudoku rules.
    """
    row, col = pos

    # Check if 'num' is already in the current row
    for j in range(9):
        if grid[row][j] == num:
            return False

    # Check if 'num' is already in the current column
    for i in range(9):
        if grid[i][col] == num:
            return False

    # Check if 'num' is already in the current 3x3 subgrid
    box_start_row = row - row % 3
    box_start_col = col - col % 3
    for i in range(3):
        for j in range(3):
            if grid[box_start_row + i][box_start_col + j] == num:
                return False

    return True

# The Sudoku puzzle grid is parsed from the visual layout in the problem description.
# A 0 represents an empty cell.
puzzle_grid = [
    [5, 0, 0, 0, 0, 8, 0, 4, 9],
    [0, 0, 0, 0, 5, 0, 0, 3, 0],
    [0, 0, 6, 3, 7, 0, 0, 0, 1],
    [1, 0, 5, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 2, 0, 8, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 1, 8],
    [7, 0, 0, 0, 0, 4, 1, 5, 0],
    [0, 0, 3, 0, 0, 2, 0, 0, 0],
    [4, 0, 9, 0, 0, 5, 0, 0, 3]
]

# Attempt to solve the puzzle
if solve_sudoku(puzzle_grid):
    # If a solution is found, extract the top horizontal line (the first row).
    top_line = puzzle_grid[0]
    
    # Print the numbers of the top line, separated by spaces.
    # The instruction "output each number in the final equation" is interpreted as
    # printing each number of the final resulting row.
    print(*top_line)
else:
    # This should not happen for a valid, solvable Sudoku puzzle.
    print("No solution exists for the given Sudoku puzzle.", file=sys.stderr)
