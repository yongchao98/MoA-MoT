def is_valid(board, row, col, num):
    # Check if the number is not in the current row, column, and 3x3 subgrid
    for i in range(9):
        if board[row][i] == num or board[i][col] == num:
            return False
    start_row, start_col = 3 * (row // 3), 3 * (col // 3)
    for i in range(start_row, start_row + 3):
        for j in range(start_col, start_col + 3):
            if board[i][j] == num:
                return False
    return True

def solve_sudoku(board):
    for row in range(9):
        for col in range(9):
            if board[row][col] == 0:
                for num in range(1, 10):
                    if is_valid(board, row, col, num):
                        board[row][col] = num
                        if solve_sudoku(board):
                            return True
                        board[row][col] = 0
                return False
    return True

# Initialize the Sudoku board with the given puzzle
sudoku_board = [
    [0, 3, 9, 0, 1, 2, 0, 0, 8],
    [8, 5, 0, 0, 0, 6, 0, 1, 0],
    [0, 0, 6, 0, 0, 0, 0, 4, 5],
    [2, 0, 3, 0, 6, 0, 0, 0, 0],
    [9, 0, 8, 2, 7, 3, 0, 0, 0],
    [0, 6, 0, 9, 0, 0, 0, 3, 0],
    [0, 7, 0, 0, 0, 0, 8, 0, 0],
    [0, 0, 0, 0, 0, 7, 0, 0, 4],
    [0, 8, 0, 4, 0, 1, 0, 0, 0]
]

# Solve the Sudoku puzzle
solve_sudoku(sudoku_board)

# Print the solved Sudoku board
for row in sudoku_board:
    print(row)