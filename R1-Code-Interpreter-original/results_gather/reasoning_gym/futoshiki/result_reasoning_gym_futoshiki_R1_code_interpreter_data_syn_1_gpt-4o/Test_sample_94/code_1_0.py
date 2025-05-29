def is_valid(board, row, col, num, horizontal_constraints, vertical_constraints):
    # Check if num is not in the current row and column
    for i in range(8):
        if board[row][i] == num or board[i][col] == num:
            return False

    # Check horizontal constraints
    if col > 0 and horizontal_constraints[row][col - 1] == '>' and board[row][col - 1] != 0 and board[row][col - 1] <= num:
        return False
    if col < 7 and horizontal_constraints[row][col] == '<' and board[row][col + 1] != 0 and board[row][col + 1] >= num:
        return False

    # Check vertical constraints
    if row > 0 and vertical_constraints[row - 1][col] == '∨' and board[row - 1][col] != 0 and board[row - 1][col] <= num:
        return False
    if row < 7 and vertical_constraints[row][col] == '∧' and board[row + 1][col] != 0 and board[row + 1][col] >= num:
        return False

    return True

def solve_futoshiki(board, horizontal_constraints, vertical_constraints):
    for row in range(8):
        for col in range(8):
            if board[row][col] == 0:
                for num in range(1, 9):
                    if is_valid(board, row, col, num, horizontal_constraints, vertical_constraints):
                        board[row][col] = num
                        if solve_futoshiki(board, horizontal_constraints, vertical_constraints):
                            return True
                        board[row][col] = 0
                return False
    return True

# Initial board setup
board = [
    [0, 1, 5, 0, 7, 0, 2, 8],
    [5, 7, 0, 1, 0, 0, 0, 3],
    [0, 0, 0, 0, 0, 8, 3, 0],
    [2, 0, 6, 0, 0, 4, 0, 0],
    [1, 0, 0, 6, 0, 0, 0, 0],
    [0, 3, 0, 0, 0, 0, 4, 1],
    [8, 2, 0, 0, 0, 5, 0, 6],
    [4, 8, 0, 0, 0, 0, 0, 5]
]

# Constraints setup
horizontal_constraints = [
    ['', '', '', '', '', '', '', ''],
    ['', '', '', '', '', '', '', ''],
    ['', '', '', '', '', '', '', ''],
    ['', '', '', '', '', '', '', ''],
    ['', '', '', '', '', '', '', ''],
    ['', '', '', '', '', '', '', ''],
    ['', '', '<', '', '', '', '', ''],
    ['', '', '', '', '', '', '', '']
]

vertical_constraints = [
    ['', '', '', '', '', '', '', ''],
    ['', '', '', '', '', '', '', ''],
    ['', '', '', '', '', '', '', ''],
    ['', '', '', '', '', '', '', '∨'],
    ['', '', '', '', '', '', '', ''],
    ['', '', '', '', '', '', '', ''],
    ['', '', '', '∨', '', '', '', ''],
    ['', '', '', '', '', '', '', '']
]

# Solve the puzzle
solve_futoshiki(board, horizontal_constraints, vertical_constraints)

# Print the solved board
for row in board:
    print(row)