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
                board[row][col] = '1'
                cols.add(col)
                diag1.add(row - col)
                diag2.add(row + col)
                if place_queens(row + 1):
                    return True
                board[row][col] = '0'
                cols.remove(col)
                diag1.remove(row - col)
                diag2.remove(row + col)
        return False

    # Initialize board and constraints
    cols = set()
    diag1 = set()
    diag2 = set()
    for r in range(8):
        for c in range(8):
            if board[r][c] == '1':
                cols.add(c)
                diag1.add(r - c)
                diag2.add(r + c)

    # Start placing queens from the first row
    place_queens(0)

    # Collect the positions of the queens
    queen_positions = []
    for r in range(8):
        for c in range(8):
            if board[r][c] == '1':
                queen_positions.append(f"{r} {c}")

    return f"<<<{', '.join(queen_positions)}>>>"

# Initial board setup
board = [
    ['0', '1', '0', '0', '0', '0', '0', '0'],
    ['0', '0', '0', '0', '0', '0', '0', '0'],
    ['0', '0', '0', '0', '0', '0', 'X', '0'],
    ['0', '0', '0', '0', '0', '0', '0', '0'],
    ['0', '0', '0', '0', '0', '0', '0', '0'],
    ['0', '0', '0', '0', '0', '0', '0', '0'],
    ['0', '0', 'X', 'X', '0', '0', '0', '0'],
    ['0', '0', '0', '0', '0', '0', '0', '0']
]

# Solve the problem and print the result
print(solve_n_queens(board))