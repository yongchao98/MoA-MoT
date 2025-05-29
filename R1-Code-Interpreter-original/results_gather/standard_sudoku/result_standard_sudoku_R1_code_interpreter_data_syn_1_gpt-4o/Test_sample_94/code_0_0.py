def is_valid(board, row, col, num):
    # Check if the number is not in the current row and column
    for i in range(9):
        if board[row][i] == num or board[i][col] == num:
            return False
    
    # Check if the number is not in the current 3x3 subgrid
    start_row, start_col = 3 * (row // 3), 3 * (col // 3)
    for i in range(3):
        for j in range(3):
            if board[start_row + i][start_col + j] == num:
                return False
    
    return True

def solve_sudoku(board):
    for row in range(9):
        for col in range(9):
            if board[row][col] == 0:  # Find an empty cell
                for num in range(1, 10):  # Try numbers 1 to 9
                    if is_valid(board, row, col, num):
                        board[row][col] = num  # Place the number
                        if solve_sudoku(board):  # Recursively solve the rest
                            return True
                        board[row][col] = 0  # Backtrack if not solvable
                return False  # No valid number found, trigger backtracking
    return True  # Solved

# Initial Sudoku grid with empty cells as 0
sudoku_grid = [
    [9, 8, 1, 0, 0, 0, 0, 7, 0],
    [0, 0, 0, 1, 0, 0, 4, 0, 0],
    [5, 0, 7, 0, 0, 0, 0, 0, 0],
    [0, 0, 9, 0, 0, 3, 0, 0, 6],
    [3, 6, 0, 0, 7, 0, 9, 0, 2],
    [0, 0, 0, 0, 6, 1, 3, 0, 0],
    [0, 0, 2, 8, 0, 0, 0, 0, 0],
    [0, 9, 5, 0, 1, 6, 0, 0, 0],
    [0, 0, 3, 0, 4, 0, 0, 0, 1]
]

solve_sudoku(sudoku_grid)
print(sudoku_grid)