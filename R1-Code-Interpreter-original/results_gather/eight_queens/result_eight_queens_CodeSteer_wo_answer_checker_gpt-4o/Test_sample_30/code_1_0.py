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
        # Skip rows that already have a queen
        if any(q[0] == r for q in queens):
            return place_queens(r + 1)
        for c in range(n):
            if can_place(r, c):
                queens.append((r, c))
                cols.add(c)
                pos_diags.add(r + c)
                neg_diags.add(r - c)
                if place_queens(r + 1):
                    return True
                # Backtrack
                queens.pop()
                cols.remove(c)
                pos_diags.remove(r + c)
                neg_diags.remove(r - c)
        return False

    # Start placing queens from the first row that doesn't have a queen
    place_queens(0)

    # Format the output
    result = ', '.join(f'{r} {c}' for r, c in sorted(queens))
    return f'<<<{result}>>>'

# Define the board
board = [
    "00000000",
    "00000000",
    "00000010",
    "00000000",
    "00000000",
    "00000000",
    "000X0000",
    "100000X0"
]

# Solve the problem
print(solve_n_queens(board))