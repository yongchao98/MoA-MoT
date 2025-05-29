def is_valid(board, row, col, num):
    # Check if num is not in the current row and column
    for i in range(7):
        if board[row][i] == num or board[i][col] == num:
            return False
    
    # Check horizontal constraints
    if col > 0:
        if (row, col - 1, '<') in constraints and board[row][col - 1] != 0 and board[row][col - 1] >= num:
            return False
        if (row, col - 1, '>') in constraints and board[row][col - 1] != 0 and board[row][col - 1] <= num:
            return False
    if col < 6:
        if (row, col, '<') in constraints and board[row][col + 1] != 0 and num >= board[row][col + 1]:
            return False
        if (row, col, '>') in constraints and board[row][col + 1] != 0 and num <= board[row][col + 1]:
            return False
    
    # Check vertical constraints
    if row > 0:
        if (row - 1, col, '∧') in constraints and board[row - 1][col] != 0 and board[row - 1][col] >= num:
            return False
        if (row - 1, col, '∨') in constraints and board[row - 1][col] != 0 and board[row - 1][col] <= num:
            return False
    if row < 6:
        if (row, col, '∧') in constraints and board[row + 1][col] != 0 and num >= board[row + 1][col]:
            return False
        if (row, col, '∨') in constraints and board[row + 1][col] != 0 and num <= board[row + 1][col]:
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
    [0, 0, 0, 4, 0, 0, 3],
    [6, 7, 4, 0, 0, 0, 0],
    [0, 4, 0, 0, 1, 3, 0],
    [0, 0, 0, 3, 0, 5, 0],
    [0, 0, 2, 0, 0, 0, 6],
    [5, 0, 0, 0, 0, 0, 0],
    [1, 6, 0, 0, 0, 0, 0]
]

# Constraints
constraints = {
    (0, 1, '<'), (0, 3, '>'), (0, 5, '>'),
    (2, 1, '<'), (2, 3, '∧'), (2, 4, '∧'),
    (3, 2, '∧'), (3, 3, '<'), (3, 4, '∧'), (3, 5, '∨'),
    (4, 2, '∨'), (4, 3, '<'), (4, 5, '<'), (4, 6, '∨'),
    (5, 2, '<'), (5, 5, '∨'),
    (6, 6, '∧')
}

solve_futoshiki(board)

# Print the solved board
for row in board:
    print(row)