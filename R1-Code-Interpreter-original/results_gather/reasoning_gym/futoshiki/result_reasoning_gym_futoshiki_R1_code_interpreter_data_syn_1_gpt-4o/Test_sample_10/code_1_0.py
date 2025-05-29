def is_valid(board, row, col, num, constraints):
    # Check if num is not in the current row and column
    for i in range(9):
        if board[row][i] == num or board[i][col] == num:
            return False

    # Check horizontal constraints
    if col > 0 and constraints[row][col - 1] == '<' and board[row][col - 1] != 0 and board[row][col - 1] >= num:
        return False
    if col > 0 and constraints[row][col - 1] == '>' and board[row][col - 1] != 0 and board[row][col - 1] <= num:
        return False
    if col < 8 and constraints[row][col] == '<' and board[row][col + 1] != 0 and board[row][col + 1] <= num:
        return False
    if col < 8 and constraints[row][col] == '>' and board[row][col + 1] != 0 and board[row][col + 1] >= num:
        return False

    # Check vertical constraints
    if row > 0 and constraints[row - 1][col] == '∧' and board[row - 1][col] != 0 and board[row - 1][col] >= num:
        return False
    if row > 0 and constraints[row - 1][col] == '∨' and board[row - 1][col] != 0 and board[row - 1][col] <= num:
        return False
    if row < 8 and constraints[row][col] == '∧' and board[row + 1][col] != 0 and board[row + 1][col] <= num:
        return False
    if row < 8 and constraints[row][col] == '∨' and board[row + 1][col] != 0 and board[row + 1][col] >= num:
        return False

    return True

def solve_futoshiki(board, constraints):
    for row in range(9):
        for col in range(9):
            if board[row][col] == 0:
                for num in range(1, 10):
                    if is_valid(board, row, col, num, constraints):
                        board[row][col] = num
                        if solve_futoshiki(board, constraints):
                            return True
                        board[row][col] = 0
                return False
    return True

# Initial board setup
board = [
    [8, 4, 0, 9, 0, 0, 2, 5, 0],
    [0, 0, 5, 3, 2, 0, 0, 1, 0],
    [3, 0, 0, 8, 0, 0, 0, 0, 0],
    [0, 0, 6, 4, 3, 2, 0, 0, 5],
    [0, 0, 2, 7, 0, 1, 9, 6, 4],
    [7, 9, 0, 0, 1, 4, 5, 0, 0],
    [0, 0, 0, 6, 0, 5, 3, 9, 8],
    [0, 0, 0, 0, 0, 0, 4, 2, 0],
    [0, 7, 9, 0, 0, 0, 0, 8, 2]
]

# Constraints setup
constraints = [
    ['', '', '', '', '', '', '', '', ''],
    ['', '', '', '', '', '', '', '', ''],
    ['', '', '', '', '', '', '', '∧', ''],
    ['', '', '', '', '', '', '', '', ''],
    ['', '', '∧', '∨', '', '', '', '', ''],
    ['', '', '', '', '', '', '', '', ''],
    ['', '', '', '', '', '', '', '', ''],
    ['', '', '', '', '', '', '', '', ''],
    ['', '', '', '', '', '', '', '', '']
]

solve_futoshiki(board, constraints)

for row in board:
    print(row)