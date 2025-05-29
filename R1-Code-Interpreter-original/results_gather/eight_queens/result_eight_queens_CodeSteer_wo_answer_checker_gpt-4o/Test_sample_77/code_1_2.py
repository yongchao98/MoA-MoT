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
                if place_queens(row + 1):
                    return True
                board[row][col] = '0'
                cols.remove(col)
                diag1.remove(row - col)
                diag2.remove(row + col)
        return False

    # Initialize the board and sets
    cols = set()
    diag1 = set()
    diag2 = set()
    queens_positions = []

    # Mark existing queens and 'X' positions
    for r in range(8):
        for c in range(8):
            if board[r][c] == '1':
                cols.add(c)
                diag1.add(r - c)
                diag2.add(r + c)
                board[r][c] = 'Q'  # Mark the initial queen
            elif board[r][c] == 'X':
                board[r][c] = 'X'

    # Start placing queens from the first row
    for r in range(8):
        if 'Q' not in board[r]:  # Skip rows with pre-placed queens
            if not place_queens(r):
                return "No solution found"

    # Collect the positions of the queens
    for r in range(8):
        for c in range(8):
            if board[r][c] == 'Q':
                queens_positions.append(f"{r} {c}")

    return f"<<<{', '.join(queens_positions)}>>>"

# Initial board setup
board = [
    ['0', '0', '0', '0', '0', '0', '0', '0'],
    ['0', '0', '0', '1', '0', '0', '0', '0'],
    ['0', '0', '0', '0', '0', '0', '0', '0'],
    ['0', '0', '0', '0', 'X', '0', '0', '0'],
    ['0', '0', '0', '0', '0', '0', '0', '0'],
    ['0', '0', '0', '0', '0', '0', '0', '0'],
    ['0', '0', 'X', 'X', '0', '0', '0', '0'],
    ['0', '0', '0', '0', '0', '0', '0', '0']
]

# Solve the problem and print the result
print(solve_n_queens(board))