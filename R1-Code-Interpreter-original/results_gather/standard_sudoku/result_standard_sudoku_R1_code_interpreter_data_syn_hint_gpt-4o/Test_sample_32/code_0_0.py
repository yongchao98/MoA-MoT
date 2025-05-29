def is_valid(board, row, col, num):
    # Check if the number is not in the current row
    for x in range(9):
        if board[row][x] == num:
            return False

    # Check if the number is not in the current column
    for x in range(9):
        if board[x][col] == num:
            return False

    # Check if the number is not in the current 3x3 subgrid
    start_row, start_col = 3 * (row // 3), 3 * (col // 3)
    for i in range(3):
        for j in range(3):
            if board[i + start_row][j + start_col] == num:
                return False

    return True

def solve_sudoku(board):
    for row in range(9):
        for col in range(9):
            if board[row][col] == 0:  # Empty cell
                for num in range(1, 10):
                    if is_valid(board, row, col, num):
                        board[row][col] = num
                        if solve_sudoku(board):
                            return True
                        board[row][col] = 0  # Backtrack
                return False
    return True

# Initialize the Sudoku board with the given puzzle
sudoku_board = [
    [9, 0, 3, 4, 0, 0, 6, 0, 0],
    [4, 0, 0, 0, 0, 0, 0, 5, 1],
    [0, 0, 0, 0, 0, 0, 0, 0, 4],
    [0, 0, 7, 8, 0, 0, 5, 6, 0],
    [0, 0, 0, 2, 6, 0, 0, 0, 0],
    [0, 0, 6, 0, 0, 1, 4, 0, 0],
    [8, 0, 1, 0, 3, 2, 0, 0, 0],
    [0, 3, 0, 5, 0, 0, 0, 9, 6],
    [0, 0, 0, 7, 9, 0, 0, 1, 2]
]

# Solve the Sudoku puzzle
solve_sudoku(sudoku_board)

# Print the solved Sudoku board
for row in sudoku_board:
    print(row)