def is_valid(board, row, col, num):
    # Check if num is not in the current row
    if num in board[row]:
        return False
    # Check if num is not in the current column
    for r in range(6):
        if board[r][col] == num:
            return False
    return True

def solve_futoshiki(board):
    for row in range(6):
        for col in range(6):
            if board[row][col] == 0:  # Empty cell
                for num in range(1, 7):  # Try numbers 1 to 6
                    if is_valid(board, row, col, num):
                        board[row][col] = num
                        if solve_futoshiki(board):
                            return True
                        board[row][col] = 0  # Backtrack
                return False
    return True

# Initial board setup with 0 representing empty cells
board = [
    [0, 3, 6, 0, 4, 0],
    [0, 6, 2, 4, 0, 0],
    [6, 2, 0, 0, 1, 4],
    [0, 0, 0, 1, 0, 0],
    [0, 4, 0, 3, 2, 0],
    [0, 0, 0, 0, 0, 0]
]

solve_futoshiki(board)
for row in board:
    print(row)