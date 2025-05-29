def is_valid(board, row, col, num):
    # Check if num is not in the current row and column
    for i in range(8):
        if board[row][i] == num or board[i][col] == num:
            return False

    # Check horizontal constraints
    if col > 0 and constraints[row][col - 1] == '<' and board[row][col - 1] != 0:
        if board[row][col - 1] >= num:
            return False
    if col < 7 and constraints[row][col] == '<' and board[row][col + 1] != 0:
        if board[row][col + 1] <= num:
            return False
    if col > 0 and constraints[row][col - 1] == '>' and board[row][col - 1] != 0:
        if board[row][col - 1] <= num:
            return False
    if col < 7 and constraints[row][col] == '>' and board[row][col + 1] != 0:
        if board[row][col + 1] >= num:
            return False

    # Check vertical constraints
    if row > 0 and constraints[row - 1][col] == '∧' and board[row - 1][col] != 0:
        if board[row - 1][col] >= num:
            return False
    if row < 7 and constraints[row][col] == '∧' and board[row + 1][col] != 0:
        if board[row + 1][col] <= num:
            return False
    if row > 0 and constraints[row - 1][col] == '∨' and board[row - 1][col] != 0:
        if board[row - 1][col] <= num:
            return False
    if row < 7 and constraints[row][col] == '∨' and board[row + 1][col] != 0:
        if board[row + 1][col] >= num:
            return False

    return True

def solve_futoshiki(board, row=0, col=0):
    if row == 8:
        return True
    if col == 8:
        return solve_futoshiki(board, row + 1, 0)
    if board[row][col] != 0:
        return solve_futoshiki(board, row, col + 1)

    for num in range(1, 9):
        if is_valid(board, row, col, num):
            board[row][col] = num
            if solve_futoshiki(board, row, col + 1):
                return True
            board[row][col] = 0

    return False

# Initial board setup with 0 representing empty cells
board = [
    [0, 7, 0, 3, 2, 4, 6, 1],
    [0, 0, 0, 0, 0, 0, 0, 3],
    [0, 0, 5, 0, 0, 1, 0, 2],
    [0, 8, 4, 0, 6, 0, 0, 0],
    [0, 2, 3, 0, 0, 0, 7, 8],
    [0, 0, 0, 0, 0, 0, 4, 0],
    [1, 6, 2, 0, 7, 0, 0, 4],
    [8, 0, 0, 2, 0, 0, 3, 0]
]

# Constraints for the puzzle
constraints = [
    ['', '<', '', '', '', '', '', ''],
    ['<', '>', '', '', '', '', '', ''],
    ['', '', '', '', '', '', '', ''],
    ['', '', '', '', '', '', '', ''],
    ['', '', '', '>', '<', '', '<', ''],
    ['', '', '', '', '', '', '<', ''],
    ['<', '>', '', '', '', '<', '', ''],
    ['', '', '>', '', '', '', '', '']
]

solve_futoshiki(board)

# Print the solved board
for row in board:
    print(row)