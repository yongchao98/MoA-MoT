def solve_n_queens(board):
    def is_safe(row, col):
        return (col not in cols and
                (row - col) not in diag1 and
                (row + col) not in diag2 and
                board[row][col] != 'X')

    def place_queens(row):
        if row == 8:
            return True
        for col in range(8):
            if is_safe(row, col):
                board[row][col] = 'Q'
                cols.add(col)
                diag1.add(row - col)
                diag2.add(row + col)
                print(f"Placing queen at {row}, {col}")
                if place_queens(row + 1):
                    return True
                print(f"Backtracking from {row}, {col}")
                board[row][col] = 0
                cols.remove(col)
                diag1.remove(row - col)
                diag2.remove(row + col)
        return False

    # Initialize the board and constraints
    cols = set()
    diag1 = set()
    diag2 = set()
    queens = []

    # Mark existing queens and constraints
    for r in range(8):
        for c in range(8):
            if board[r][c] == 1:
                cols.add(c)
                diag1.add(r - c)
                diag2.add(r + c)
                queens.append((r, c))
                board[r][c] = 'Q'  # Mark the queen on the board

    # Solve the board
    place_queens(0)

    # Collect the positions of the queens
    for r in range(8):
        for c in range(8):
            if board[r][c] == 'Q' and (r, c) not in queens:
                queens.append((r, c))

    # Format the output
    result = ', '.join(f'{r} {c}' for r, c in sorted(queens))
    return f'<<<{result}>>>'

# Define the board
board = [
    [0, 0, 0, 0, 0, 0, 0, 0],
    [1, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 'X', 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 'X', 0, 0, 'X'],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, 0, 0, 0, 0]
]

# Solve and print the result
print(solve_n_queens(board))