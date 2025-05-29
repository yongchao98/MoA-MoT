def is_valid(board, row, col, num):
    # Check if num is not in the current row and column
    for i in range(7):
        if board[row][i] == num or board[i][col] == num:
            return False
    
    # Check horizontal constraints
    if col > 0 and board[row][col - 1] != 0:
        if (row, col - 1) == (0, 2) and board[row][col - 1] >= num:
            return False
    if col < 6 and board[row][col + 1] != 0:
        if (row, col) == (0, 2) and board[row][col + 1] <= num:
            return False
    
    # Check vertical constraints
    if row > 0 and board[row - 1][col] != 0:
        if (row - 1, col) == (4, 1) and board[row - 1][col] >= num:
            return False
    if row < 6 and board[row + 1][col] != 0:
        if (row, col) == (4, 1) and board[row + 1][col] <= num:
            return False
    
    return True

def solve_futoshiki(board):
    for row in range(7):
        for col in range(7):
            if board[row][col] == 0:
                for num in range(1, 8):
                    if is_valid(board, row, col, num):
                        board[row][col] = num
                        if solve_futoshiki(board):
                            return True
                        board[row][col] = 0
                return False
    return True

# Initial board setup
board = [
    [0, 6, 4, 0, 0, 0, 2],
    [7, 0, 0, 0, 1, 0, 0],
    [0, 2, 0, 0, 4, 6, 0],
    [4, 0, 2, 0, 0, 1, 0],
    [2, 0, 6, 0, 3, 0, 0],
    [0, 1, 3, 7, 0, 0, 0],
    [0, 0, 0, 3, 6, 0, 0]
]

solve_futoshiki(board)

# Print the solved board
for row in board:
    print(row)