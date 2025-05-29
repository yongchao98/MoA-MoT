def solve_n_queens(board):
    n = 8
    queens = []
    restricted = set()
    cols = set()
    diag1 = set()
    diag2 = set()

    # Parse the board
    for r in range(n):
        for c in range(n):
            if board[r][c] == '1':
                queens.append((r, c))
                cols.add(c)
                diag1.add(r - c)
                diag2.add(r + c)
            elif board[r][c] == 'X':
                restricted.add((r, c))

    def place_queens(row):
        if row == n:
            return True
        for col in range(n):
            if col in cols or (row - col) in diag1 or (row + col) in diag2 or (row, col) in restricted:
                continue
            # Place queen
            queens.append((row, col))
            cols.add(col)
            diag1.add(row - col)
            diag2.add(row + col)
            # Recurse
            if place_queens(row + 1):
                return True
            # Backtrack
            queens.pop()
            cols.remove(col)
            diag1.remove(row - col)
            diag2.remove(row + col)
        return False

    place_queens(len(queens))

    # Format the output
    result = ', '.join(f'{r} {c}' for r, c in sorted(queens))
    return f'<<<{result}>>>'

# Example board
board = [
    "00001000",
    "01000000",
    "00000000",
    "00000000",
    "00000000",
    "00X00000",
    "00X00000",
    "00000000"
]

print(solve_n_queens(board))