def solve_n_queens(board):
    def is_valid(row, col):
        return (col not in cols and
                (row - col) not in diag1 and
                (row + col) not in diag2 and
                board[row][col] != 'X')

    def place_queens(row):
        if row == 8:
            return True
        # Skip rows that already have a queen
        if any(board[row][c] == '1' for c in range(8)):
            return place_queens(row + 1)
        for col in range(8):
            if is_valid(row, col):
                print(f"Placing queen at ({row}, {col})")
                queens.append((row, col))
                cols.add(col)
                diag1.add(row - col)
                diag2.add(row + col)
                if place_queens(row + 1):
                    return True
                print(f"Backtracking from ({row}, {col})")
                queens.pop()
                cols.remove(col)
                diag1.remove(row - col)
                diag2.remove(row + col)
        return False

    queens = []
    cols = set()
    diag1 = set()
    diag2 = set()

    # Initialize the board and mark existing queens
    for r in range(8):
        for c in range(8):
            if board[r][c] == '1':
                print(f"Initial queen at ({r}, {c})")
                queens.append((r, c))
                cols.add(c)
                diag1.add(r - c)
                diag2.add(r + c)

    # Start placing queens from the first row
    place_queens(0)

    # Format the output
    result = ', '.join(f'{r} {c}' for r, c in sorted(queens))
    return f'<<<{result}>>>'

# Define the board
board = [
    "00000000",
    "00000000",
    "00000X01",
    "00000000",
    "00000000",
    "00000000",
    "00000000",
    "00000000"
]

# Solve the problem
print(solve_n_queens(board))