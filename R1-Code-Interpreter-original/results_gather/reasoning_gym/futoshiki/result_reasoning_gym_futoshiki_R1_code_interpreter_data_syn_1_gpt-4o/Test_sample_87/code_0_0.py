def is_valid(board, row, col, num):
    # Check if num is not in the current row and column
    for i in range(6):
        if board[row][i] == num or board[i][col] == num:
            return False

    # Check inequality constraints
    if row == 0 and col == 1 and not (num > board[1][1]):
        return False
    if row == 1 and col == 5 and not (num < board[2][5]):
        return False
    if row == 2 and col == 2 and not (num < board[2][3]):
        return False
    if row == 3 and col == 4 and not (num < board[3][5]):
        return False
    if row == 5 and col == 1 and not (num > board[5][2]):
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
    [2, 3, 0, 0, 0, 6],
    [0, 6, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 2],
    [0, 2, 0, 0, 3, 5],
    [0, 0, 5, 0, 0, 1],
    [4, 0, 0, 0, 1, 0]
]

solve_futoshiki(board)

# Print the solved board
for row in board:
    print(row)