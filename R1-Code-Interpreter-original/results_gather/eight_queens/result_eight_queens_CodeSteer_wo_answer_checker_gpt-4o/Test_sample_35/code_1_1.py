def solve_queens(board):
    n = 8
    queens = []
    occupied_cols = set()
    occupied_diag1 = set()  # r - c
    occupied_diag2 = set()  # r + c

    # Mark the already placed queens and blocked positions
    for r in range(n):
        for c in range(n):
            if board[r][c] == '1':
                queens.append((r, c))
                occupied_cols.add(c)
                occupied_diag1.add(r - c)
                occupied_diag2.add(r + c)
            elif board[r][c] == 'X':
                occupied_cols.add(c)
                occupied_diag1.add(r - c)
                occupied_diag2.add(r + c)

    def backtrack(row):
        if row == n:
            return True
        for col in range(n):
            if (col not in occupied_cols and
                (row - col) not in occupied_diag1 and
                (row + col) not in occupied_diag2 and
                board[row][col] != 'X'):
                # Place the queen
                queens.append((row, col))
                occupied_cols.add(col)
                occupied_diag1.add(row - col)
                occupied_diag2.add(row + col)

                # Recurse to place the next queen
                if backtrack(row + 1):
                    return True

                # Backtrack
                queens.pop()
                occupied_cols.remove(col)
                occupied_diag1.remove(row - col)
                occupied_diag2.remove(row + col)
        return False

    # Start the backtracking from the first row
    backtrack(0)

    # Format the output
    result = ', '.join(f'{r} {c}' for r, c in sorted(queens))
    return f'<<<{result}>>>'

# Define the board
board = [
    "00000000",
    "10000000",
    "00010000",
    "0000X000",
    "00000000",
    "00000000",
    "0000000X",
    "000X0000"
]

# Solve and print the result
print(solve_queens(board))