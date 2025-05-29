def is_safe(board, row, col, n):
    # Check this column on upper side
    for i in range(row):
        if board[i][col] == 1:
            return False

    # Check upper diagonal on left side
    for i, j in zip(range(row, -1, -1), range(col, -1, -1)):
        if board[i][j] == 1:
            return False

    # Check upper diagonal on right side
    for i, j in zip(range(row, -1, -1), range(col, n)):
        if board[i][j] == 1:
            return False

    return True

def solve_n_queens(board, row, n):
    if row == n:
        return True

    for col in range(n):
        if board[row][col] == 'X' or board[row][col] == 1:
            continue

        if is_safe(board, row, col, n):
            board[row][col] = 1
            print(f"Placing queen at {row}, {col}")
            print_board(board)

            if solve_n_queens(board, row + 1, n):
                return True

            # Backtrack
            board[row][col] = 0
            print(f"Removing queen from {row}, {col}")
            print_board(board)

    return False

def print_board(board):
    for row in board:
        print(' '.join(str(x) for x in row))
    print()

def find_queen_positions(board):
    n = len(board)
    if solve_n_queens(board, 0, n):
        positions = []
        for i in range(n):
            for j in range(n):
                if board[i][j] == 1:
                    positions.append(f"{i} {j}")
        return positions
    else:
        return "No solution found"

# Initial board setup
board = [
    [0, 0, 0, 0, 0, 0, 0, 0],
    [1, 0, 0, 0, 0, 0, 0, 0],
    ['X', 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0]
]

positions = find_queen_positions(board)
print(f"<<<{', '.join(positions)}>>>")