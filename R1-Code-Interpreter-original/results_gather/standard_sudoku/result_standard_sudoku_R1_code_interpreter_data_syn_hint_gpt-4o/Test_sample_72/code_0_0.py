def is_valid(board, row, col, num):
    # Check if the number is not in the current row, column, and 3x3 subgrid
    for i in range(9):
        if board[row][i] == num or board[i][col] == num:
            return False
    start_row, start_col = 3 * (row // 3), 3 * (col // 3)
    for i in range(3):
        for j in range(3):
            if board[start_row + i][start_col + j] == num:
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

# Initial Sudoku grid with 0 representing empty cells
sudoku_grid = [
    [0, 0, 0, 1, 3, 5, 0, 7, 9],
    [0, 0, 0, 6, 9, 0, 1, 0, 3],
    [0, 0, 0, 0, 0, 0, 5, 0, 6],
    [2, 3, 0, 0, 0, 0, 7, 0, 0],
    [7, 0, 0, 0, 1, 4, 0, 0, 2],
    [0, 0, 0, 0, 0, 0, 0, 4, 8],
    [0, 0, 6, 8, 0, 0, 0, 0, 0],
    [5, 0, 1, 0, 2, 6, 0, 0, 0],
    [0, 2, 8, 0, 0, 0, 9, 0, 4]
]

solve_sudoku(sudoku_grid)
print(sudoku_grid)