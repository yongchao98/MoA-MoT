def is_valid(board, row, col, num):
    # Check if num is not in the current row and column
    for i in range(7):
        if board[row][i] == num or board[i][col] == num:
            return False

    # Check horizontal constraints
    if col > 0 and board[row][col - 1] != 0:
        if (board[row][col - 1] < num and constraints[row][col - 1] == '>') or \
           (board[row][col - 1] > num and constraints[row][col - 1] == '<'):
            return False
    if col < 6 and board[row][col + 1] != 0:
        if (board[row][col + 1] < num and constraints[row][col] == '<') or \
           (board[row][col + 1] > num and constraints[row][col] == '>'):
            return False

    # Check vertical constraints
    if row > 0 and board[row - 1][col] != 0:
        if (board[row - 1][col] < num and constraints[row - 1][col] == '∨') or \
           (board[row - 1][col] > num and constraints[row - 1][col] == '∧'):
            return False
    if row < 6 and board[row + 1][col] != 0:
        if (board[row + 1][col] < num and constraints[row][col] == '∧') or \
           (board[row + 1][col] > num and constraints[row][col] == '∨'):
            return False

    return True

def solve_futoshiki(board, row=0, col=0):
    if row == 7:
        return True
    if col == 7:
        return solve_futoshiki(board, row + 1, 0)
    if board[row][col] != 0:
        return solve_futoshiki(board, row, col + 1)

    for num in range(1, 8):
        if is_valid(board, row, col, num):
            board[row][col] = num
            if solve_futoshiki(board, row, col + 1):
                return True
            board[row][col] = 0

    return False

# Initial board setup
board = [
    [3, 0, 0, 0, 4, 0, 0],
    [4, 0, 1, 0, 0, 6, 0],
    [7, 0, 0, 0, 5, 0, 0],
    [6, 4, 0, 3, 7, 2, 0],
    [0, 0, 0, 0, 3, 4, 0],
    [0, 1, 0, 7, 0, 0, 0],
    [0, 5, 6, 0, 0, 7, 0]
]

# Constraints setup
constraints = [
    ['', '', '', '', '', '', ''],
    ['', '', '', '', '', '', ''],
    ['', '', '', '', '', '', ''],
    ['', '', '', '', '', '', ''],
    ['', '', '', '>', '', '', ''],
    ['', '', '', '', '', '', ''],
    ['', '', '', '', '', '', '']
]

solve_futoshiki(board)

for row in board:
    print(row)