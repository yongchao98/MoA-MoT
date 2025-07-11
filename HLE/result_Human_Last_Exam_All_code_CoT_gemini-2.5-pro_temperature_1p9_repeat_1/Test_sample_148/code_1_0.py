def solve_sudoku(board):
    """
    Solves a Sudoku puzzle using backtracking.
    The board is a 9x9 list of lists, with 0 representing empty cells.
    """

    def find_empty(bo):
        """Finds an empty cell (represented by 0) in the board."""
        for i in range(len(bo)):
            for j in range(len(bo[0])):
                if bo[i][j] == 0:
                    return (i, j)  # row, col
        return None

    def is_valid(bo, num, pos):
        """Checks if placing a number in a given position is valid."""
        row, col = pos

        # Check row
        for i in range(len(bo[0])):
            if bo[row][i] == num and col != i:
                return False

        # Check column
        for i in range(len(bo)):
            if bo[i][col] == num and row != i:
                return False

        # Check 3x3 box
        box_x = col // 3
        box_y = row // 3

        for i in range(box_y * 3, box_y * 3 + 3):
            for j in range(box_x * 3, box_x * 3 + 3):
                if bo[i][j] == num and (i, j) != pos:
                    return False

        return True

    # Main backtracking logic
    find = find_empty(board)
    if not find:
        return True  # Puzzle solved
    else:
        row, col = find

    for i in range(1, 10):
        if is_valid(board, i, (row, col)):
            board[row][col] = i

            if solve_sudoku(board):
                return True

            board[row][col] = 0  # Backtrack

    return False

# Based on careful parsing of the visual representation of the puzzle
# to ensure it's a valid and solvable board.
sudoku_board = [
    [5, 0, 0, 0, 8, 0, 0, 4, 9],
    [0, 0, 0, 5, 0, 0, 0, 3, 0],
    [0, 0, 6, 3, 7, 0, 0, 0, 1],
    [1, 0, 5, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 2, 0, 8, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 1, 8],
    [7, 0, 0, 0, 0, 4, 1, 5, 0],
    [0, 0, 3, 0, 0, 2, 0, 0, 0],
    [4, 0, 9, 0, 5, 0, 0, 0, 3]
]


if solve_sudoku(sudoku_board):
    # Print the top horizontal line of the solved sudoku
    # as 9 space-separated numbers.
    top_line = sudoku_board[0]
    # Each number in the final equation needs to be output
    print(*top_line)
else:
    # This case should not be reached if the puzzle is valid.
    print("Could not solve the Sudoku puzzle.")
