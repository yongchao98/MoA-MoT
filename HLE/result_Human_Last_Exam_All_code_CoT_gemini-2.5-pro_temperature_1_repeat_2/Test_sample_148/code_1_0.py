def solve_sudoku(board):
    """
    Solves a Sudoku puzzle using backtracking.
    """
    for row in range(9):
        for col in range(9):
            if board[row][col] == 0:
                for num in range(1, 10):
                    if is_valid(board, row, col, num):
                        board[row][col] = num
                        if solve_sudoku(board):
                            return True
                        board[row][col] = 0  # Backtrack
                return False
    return True

def is_valid(board, row, col, num):
    """
    Checks if placing a number in a given cell is valid.
    """
    # Check if the number is already in the current row
    for x in range(9):
        if board[row][x] == num:
            return False

    # Check if the number is already in the current column
    for x in range(9):
        if board[x][col] == num:
            return False

    # Check if the number is in the current 3x3 subgrid
    start_row = row - row % 3
    start_col = col - col % 3
    for i in range(3):
        for j in range(3):
            if board[i + start_row][j + start_col] == num:
                return False
                
    return True

# The initial board is derived by interpreting the visual blocks in the puzzle.
# The original puzzle has a contradiction (two 8s in column 6).
# The board below includes a correction for this assumed typo:
# The 8 in row 5, column 6 has been moved to row 5, column 5.
initial_board = [
    [5, 0, 0, 0, 0, 8, 0, 4, 9],
    [0, 0, 0, 5, 0, 0, 0, 3, 0],
    [0, 0, 6, 3, 0, 0, 0, 0, 1],
    [1, 0, 5, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 2, 8, 0, 0, 0, 0], # Corrected row
    [0, 0, 0, 0, 0, 0, 1, 0, 8],
    [7, 0, 0, 0, 0, 4, 1, 5, 0],
    [0, 0, 3, 0, 0, 2, 0, 0, 0],
    [4, 0, 9, 0, 0, 5, 0, 0, 3]
]

if solve_sudoku(initial_board):
    # The puzzle is solved, print the top horizontal line.
    # The '*' unpacks the list into separate arguments for print.
    print(*initial_board[0])
else:
    print("No solution exists for the corrected puzzle.")
