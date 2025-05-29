def solve_n_queens(board):
    def is_safe(row, col):
        return not (cols[col] or diag1[row - col] or diag2[row + col] or board[row][col] == 'X')

    def place_queen(row):
        if row == 8:
            return True
        for col in range(8):
            if is_safe(row, col):
                queens.append((row, col))
                cols[col] = diag1[row - col] = diag2[row + col] = True
                if place_queen(row + 1):
                    return True
                queens.pop()
                cols[col] = diag1[row - col] = diag2[row + col] = False
        return False

    # Initialize board state
    queens = []
    cols = [False] * 8
    diag1 = [False] * 15
    diag2 = [False] * 15

    # Mark existing queens and 'X' positions
    for r in range(8):
        for c in range(8):
            if board[r][c] == 1:
                queens.append((r, c))
                cols[c] = diag1[r - c] = diag2[r + c] = True

    # Start placing queens from the first row
    place_queen(0)

    # Format the output
    result = ', '.join(f'{r} {c}' for r, c in queens)
    print(f'<<<{result}>>>')

# Initial board setup
board = [
    [0, 0, 0, 0, 0, 0, 0, 0],
    [1, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 1, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0]
]

solve_n_queens(board)