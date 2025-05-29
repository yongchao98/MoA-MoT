def is_valid(board, row, col, num, constraints):
    # Check if num is not in the current row and column
    for i in range(7):
        if board[row][i] == num or board[i][col] == num:
            return False

    # Check inequality constraints
    for (r1, c1, r2, c2) in constraints:
        if (r1, c1) == (row, col) and board[r2][c2] != 0 and not (num < board[r2][c2]):
            return False
        if (r2, c2) == (row, col) and board[r1][c1] != 0 and not (board[r1][c1] < num):
            return False

    return True

def solve_futoshiki(board, constraints):
    for row in range(7):
        for col in range(7):
            if board[row][col] == 0:
                for num in range(1, 8):
                    if is_valid(board, row, col, num, constraints):
                        board[row][col] = num
                        if solve_futoshiki(board, constraints):
                            return True
                        board[row][col] = 0
                return False
    return True

# Initial board setup with constraints
board = [
    [0, 0, 5, 0, 0, 0, 0],
    [0, 0, 0, 4, 1, 0, 0],
    [0, 1, 0, 0, 0, 7, 2],
    [0, 0, 1, 6, 2, 0, 4],
    [0, 0, 0, 0, 0, 0, 1],
    [0, 0, 0, 0, 5, 0, 0],
    [0, 6, 0, 0, 7, 0, 0]
]

# Constraints: (row1, col1, row2, col2) means board[row1][col1] < board[row2][col2]
constraints = [
    (0, 5, 0, 6), (1, 2, 1, 3), (2, 2, 3, 2), (3, 2, 3, 3),
    (5, 3, 5, 4), (6, 0, 6, 1), (6, 5, 6, 6), (6, 1, 5, 1),
    (6, 4, 6, 5), (1, 3, 2, 3), (3, 3, 4, 3)
]

solve_futoshiki(board, constraints)

# Print the solved board
for row in board:
    print(row)