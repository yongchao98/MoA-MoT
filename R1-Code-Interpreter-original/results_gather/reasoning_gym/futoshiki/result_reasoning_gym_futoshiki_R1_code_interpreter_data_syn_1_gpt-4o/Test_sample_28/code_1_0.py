def is_valid(board, row, col, num):
    # Check if num is not in the current row and column
    for i in range(6):
        if board[row][i] == num or board[i][col] == num:
            return False
    
    # Check horizontal constraints
    if col > 0 and constraints[row][col - 1] == '<' and board[row][col - 1] != 0:
        if not (board[row][col - 1] < num):
            return False
    if col > 0 and constraints[row][col - 1] == '>' and board[row][col - 1] != 0:
        if not (board[row][col - 1] > num):
            return False
    if col < 5 and constraints[row][col] == '<' and board[row][col + 1] != 0:
        if not (num < board[row][col + 1]):
            return False
    if col < 5 and constraints[row][col] == '>' and board[row][col + 1] != 0:
        if not (num > board[row][col + 1]):
            return False
    
    # Check vertical constraints
    if row > 0 and constraints[row - 1][col] == '∧' and board[row - 1][col] != 0:
        if not (board[row - 1][col] < num):
            return False
    if row > 0 and constraints[row - 1][col] == '∨' and board[row - 1][col] != 0:
        if not (board[row - 1][col] > num):
            return False
    if row < 5 and constraints[row][col] == '∧' and board[row + 1][col] != 0:
        if not (num < board[row + 1][col]):
            return False
    if row < 5 and constraints[row][col] == '∨' and board[row + 1][col] != 0:
        if not (num > board[row + 1][col]):
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

# Initial board setup with 0 representing blanks
board = [
    [4, 0, 1, 0, 0, 0],
    [0, 0, 5, 0, 0, 0],
    [0, 0, 0, 4, 1, 0],
    [0, 0, 0, 0, 5, 4],
    [1, 0, 0, 6, 0, 3],
    [0, 0, 3, 0, 0, 0]
]

# Constraints for each cell
constraints = [
    ['', '', '<', '', '', ''],
    ['<', '', '', '<', '', '<', ''],
    ['', '', '', '>', '', ''],
    ['', '', '<', '', '', ''],
    ['<', '', '', '', '', ''],
    ['', '', '', '<', '', '']
]

solve_futoshiki(board)
for row in board:
    print(row)