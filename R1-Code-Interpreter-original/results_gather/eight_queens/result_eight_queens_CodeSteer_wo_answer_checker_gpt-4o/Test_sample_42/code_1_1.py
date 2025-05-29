def solve_n_queens(board):
    def is_safe(row, col):
        return not (cols[col] or diag1[row - col] or diag2[row + col])

    def place_queen(row):
        if row == 8:
            return True
        for col in range(8):
            if board[row][col] == 'X' or not is_safe(row, col):
                continue
            # Place the queen
            queens.append((row, col))
            cols[col] = diag1[row - col] = diag2[row + col] = True
            print(f"Placing queen at: {row}, {col}")
            if place_queen(row + 1):
                return True
            # Backtrack
            print(f"Removing queen from: {row}, {col}")
            queens.pop()
            cols[col] = diag1[row - col] = diag2[row + col] = False
        return False

    # Initialize the board
    queens = []
    cols = [False] * 8
    diag1 = [False] * 15  # 8 + 7 (row - col can range from -7 to 7)
    diag2 = [False] * 15  # 8 + 7 (row + col can range from 0 to 14)

    # Mark the existing queens
    for r in range(8):
        for c in range(8):
            if board[r][c] == '1':
                queens.append((r, c))
                cols[c] = diag1[r - c] = diag2[r + c] = True

    # Start placing queens from the first row
    place_queen(0)

    # Format the output
    result = ', '.join(f'{r} {c}' for r, c in queens)
    return f'<<<{result}>>>'

# Example board
board = [
    "01000000",
    "00000000",
    "00000000",
    "00000010",
    "00000000",
    "00000001",
    "00000000",
    "00000000"
]

print(solve_n_queens(board))