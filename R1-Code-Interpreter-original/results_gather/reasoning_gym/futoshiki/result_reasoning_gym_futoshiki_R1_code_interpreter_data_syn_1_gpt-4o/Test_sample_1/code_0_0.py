def is_valid(board, row, col, num):
    # Check if num is not in the current row and column
    for i in range(8):
        if board[row][i] == num or board[i][col] == num:
            return False
    
    # Check inequality constraints
    if row > 0 and board[row-1][col] != 0:
        if constraints[row-1][col] == '∧' and not (board[row-1][col] < num):
            return False
        if constraints[row-1][col] == '∨' and not (board[row-1][col] > num):
            return False
    if col > 0 and board[row][col-1] != 0:
        if constraints[row][col-1] == '<' and not (board[row][col-1] < num):
            return False
        if constraints[row][col-1] == '>' and not (board[row][col-1] > num):
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

# Initial board setup
board = [
    [1, 5, 0, 0, 0, 3, 0, 0],
    [0, 0, 2, 6, 0, 0, 3, 0],
    [7, 0, 4, 0, 0, 0, 8, 1],
    [2, 8, 5, 0, 0, 0, 0, 7],
    [0, 0, 0, 0, 3, 0, 1, 6],
    [0, 4, 0, 0, 0, 0, 0, 5],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 8, 4, 5, 7, 0, 0]
]

# Constraints setup
constraints = [
    ['', '', '', '', '', '', '', ''],
    ['', '', '', '', '', '', '', ''],
    ['', '', '', '', '', '', '', ''],
    ['', '', '', '', '', '', '', ''],
    ['', '', '', '', '<', '>', '', ''],
    ['∧', '', '', '', '', '', '', ''],
    ['∨', '', '', '', '', '', '', ''],
    ['', '', '', '', '', '', '', '']
]

solve_futoshiki(board)

# Print the solved board
for row in board:
    print(row)