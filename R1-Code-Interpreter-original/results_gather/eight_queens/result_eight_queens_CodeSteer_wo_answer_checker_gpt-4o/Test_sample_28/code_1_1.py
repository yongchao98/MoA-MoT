def solve_n_queens(board):
    n = 8
    queens = []
    cols = set()
    pos_diags = set()  # r + c
    neg_diags = set()  # r - c

    # Parse the board to find existing queens and restricted positions
    for r in range(n):
        for c in range(n):
            if board[r][c] == '1':
                queens.append((r, c))
                cols.add(c)
                pos_diags.add(r + c)
                neg_diags.add(r - c)

    def can_place(r, c):
        return c not in cols and (r + c) not in pos_diags and (r - c) not in neg_diags and board[r][c] != 'X'

    def place_queens(r):
        if r == n:
            return True
        for c in range(n):
            if can_place(r, c):
                print(f"Placing queen at ({r}, {c})")
                queens.append((r, c))
                cols.add(c)
                pos_diags.add(r + c)
                neg_diags.add(r - c)
                if place_queens(r + 1):
                    return True
                # Backtrack
                print(f"Removing queen from ({r}, {c})")
                queens.pop()
                cols.remove(c)
                pos_diags.remove(r + c)
                neg_diags.remove(r - c)
        return False

    # Start placing queens from the first row
    place_queens(0)

    # Format the output
    result = ', '.join(f'{r} {c}' for r, c in sorted(queens))
    return f'<<<{result}>>>'

# Define the board
board = [
    "00X00000",
    "01000000",
    "00000000",
    "0000000X",
    "00000000",
    "00000000",
    "00001000",
    "00000000"
]

# Solve the problem
print(solve_n_queens(board))