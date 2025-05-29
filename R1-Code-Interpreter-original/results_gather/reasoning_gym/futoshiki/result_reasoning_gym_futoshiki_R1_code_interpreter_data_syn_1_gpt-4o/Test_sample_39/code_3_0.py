def is_valid(board, row, col, num):
    # Check if num is not in the current row and column
    for i in range(6):
        if board[row][i] == num or board[i][col] == num:
            return False
    
    # Check horizontal constraints
    if col < 5 and board[row][col + 1] != 0 and board[row][col + 1] < num:
        return False
    if col > 0 and board[row][col - 1] != 0 and board[row][col - 1] > num:
        return False
    
    # Check vertical constraints
    if row < 5 and board[row + 1][col] != 0 and board[row + 1][col] < num:
        return False
    if row > 0 and board[row - 1][col] != 0 and board[row - 1][col] > num:
        return False
    
    return True

def solve_futoshiki(board):
    for row in range(6):
        for col in range(6):
            if board[row][col] == 0:
                for num in range(1, 7):
                    if is_valid(board, row, col, num):
                        board[row][col] = num
                        if solve_futoshiki(board):
                            return True
                        board[row][col] = 0
                return False
    return True

# Initial board setup
board = [
    [6, 0, 4, 0, 2, 0],
    [0, 6, 0, 0, 5, 2],
    [0, 0, 6, 3, 0, 0],
    [5, 2, 0, 0, 3, 0],
    [0, 0, 0, 4, 0, 1],
    [0, 0, 5, 0, 0, 0]
]

solve_futoshiki(board)
for row in board:
    print(row)