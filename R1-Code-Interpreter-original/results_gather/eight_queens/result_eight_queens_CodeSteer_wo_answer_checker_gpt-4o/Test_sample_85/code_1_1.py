def solve_queens(board):
    n = 8
    queens = []
    solutions = []
    occupied_cols = set()
    occupied_diag1 = set()  # r - c
    occupied_diag2 = set()  # r + c

    def is_safe(r, c):
        return (c not in occupied_cols and
                (r - c) not in occupied_diag1 and
                (r + c) not in occupied_diag2 and
                board[r][c] != 'X')

    def place_queens(row):
        if row == n:
            solutions.append(queens.copy())
            return True

        for col in range(n):
            if is_safe(row, col):
                queens.append((row, col))
                occupied_cols.add(col)
                occupied_diag1.add(row - col)
                occupied_diag2.add(row + col)

                if place_queens(row + 1):
                    return True

                # Backtrack
                queens.pop()
                occupied_cols.remove(col)
                occupied_diag1.remove(row - col)
                occupied_diag2.remove(row + col)

        return False

    # Start placing queens from the first row
    place_queens(0)

    # Format the output
    if solutions:
        result = ', '.join(f'{r} {c}' for r, c in sorted(solutions[0]))
        return f'<<<{result}>>>'
    else:
        return "No solution found"

# Define the board
board = [
    "0X100000",
    "00000000",
    "00000000",
    "0X000000",
    "00000000",
    "000000X0",
    "00000000",
    "00000000"
]

# Solve and print the result
print(solve_queens(board))