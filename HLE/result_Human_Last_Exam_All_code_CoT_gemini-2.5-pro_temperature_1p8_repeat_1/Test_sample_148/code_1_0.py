import sys

def solve_sudoku(board):
    """
    Solves a Sudoku puzzle using backtracking.
    The board is modified in-place.
    """
    pos = find_empty_cell(board)
    if not pos:
        return True  # Puzzle solved

    row, col = pos
    for num in range(1, 10):
        if is_safe(board, row, col, num):
            board[row][col] = num
            if solve_sudoku(board):
                return True
            board[row][col] = 0  # Backtrack
    return False

def find_empty_cell(board):
    """Finds an empty cell (represented by 0) in the board."""
    for i in range(9):
        for j in range(9):
            if board[i][j] == 0:
                return (i, j)
    return None

def is_safe(board, row, col, num):
    """Checks if it's safe to place a number in a given cell."""
    # Check row for the same number
    for x in range(9):
        if board[row][x] == num:
            return False
    # Check column for the same number
    for x in range(9):
        if board[x][col] == num:
            return False
    # Check 3x3 box for the same number
    start_row = row - row % 3
    start_col = col - col % 3
    for i in range(3):
        for j in range(3):
            if board[i + start_row][j + start_col] == num:
                return False
    return True

def main():
    """
    Main function to define, solve, and print the Sudoku result.
    """
    # The Sudoku grid parsed from the visual representation in the prompt.
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

    if solve_sudoku(grid):
        # Print the top horizontal line of the solved sudoku.
        # This will be the first row of the grid.
        # The elements are converted to strings and joined with spaces.
        print(' '.join(map(str, grid[0])))
    else:
        print("No solution exists for the given Sudoku puzzle.")

if __name__ == "__main__":
    main()
