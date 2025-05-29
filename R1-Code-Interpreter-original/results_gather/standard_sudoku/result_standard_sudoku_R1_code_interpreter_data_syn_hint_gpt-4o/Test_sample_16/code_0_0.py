def is_valid(board, row, col, num):
    # Check if num is not in the current row
    for x in range(9):
        if board[row][x] == num:
            return False

    # Check if num is not in the current column
    for x in range(9):
        if board[x][col] == num:
            return False

    # Check if num is not in the current 3x3 subgrid
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

# Initial Sudoku grid
sudoku_grid = [
    [0, 0, 0, 4, 0, 6, 7, 0, 0],
    [1, 0, 4, 0, 9, 0, 0, 5, 0],
    [9, 7, 0, 3, 0, 8, 0, 0, 0],
    [2, 4, 8, 0, 0, 0, 0, 0, 0],
    [0, 5, 9, 0, 0, 0, 0, 0, 0],
    [7, 0, 0, 0, 0, 9, 4, 0, 2],
    [0, 3, 0, 0, 0, 0, 0, 0, 9],
    [0, 6, 0, 0, 0, 1, 0, 0, 7],
    [0, 9, 1, 2, 0, 3, 0, 8, 0]
]

solve_sudoku(sudoku_grid)
print(sudoku_grid)