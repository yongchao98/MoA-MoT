def is_safe(board, row, col):
    # Check if the column is safe
    for r in range(row):
        if board[r][col] == 1:
            return False

    # Check the upper left diagonal
    for r, c in zip(range(row, -1, -1), range(col, -1, -1)):
        if board[r][c] == 1:
            return False

    # Check the upper right diagonal
    for r, c in zip(range(row, -1, -1), range(col, len(board))):
        if board[r][c] == 1:
            return False

    return True

def solve_n_queens(board, row, queens):
    if row == len(board):
        return True

    for col in range(len(board)):
        if board[row][col] == 0 and is_safe(board, row, col):
            board[row][col] = 1
            queens.append((row, col))
            if solve_n_queens(board, row + 1, queens):
                return True
            board[row][col] = 0
            queens.pop()

    return False

def place_queens_on_board():
    board = [
        [0, 0, 0, 0, 0, 0, -1, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [-1, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, -1]
    ]

    queens = []
    solve_n_queens(board, 0, queens)

    # Format the output
    result = ', '.join(f"{r} {c}" for r, c in queens)
    print(f"<<<{result}>>>")

place_queens_on_board()