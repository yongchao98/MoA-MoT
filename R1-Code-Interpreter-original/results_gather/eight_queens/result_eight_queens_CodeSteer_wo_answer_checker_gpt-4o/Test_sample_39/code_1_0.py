def solve_n_queens(board):
    n = 8
    queens = []
    columns = set()
    pos_diagonals = set()  # r + c
    neg_diagonals = set()  # r - c

    # Parse the board to find existing queens and blocked positions
    for r in range(n):
        for c in range(n):
            if board[r][c] == '1':
                queens.append((r, c))
                columns.add(c)
                pos_diagonals.add(r + c)
                neg_diagonals.add(r - c)

    def is_valid(r, c):
        return (c not in columns and
                (r + c) not in pos_diagonals and
                (r - c) not in neg_diagonals and
                board[r][c] != 'X')

    def place_queens(r):
        if r == n:
            return True
        for c in range(n):
            if is_valid(r, c):
                queens.append((r, c))
                columns.add(c)
                pos_diagonals.add(r + c)
                neg_diagonals.add(r - c)

                if place_queens(r + 1):
                    return True

                # Backtrack
                queens.pop()
                columns.remove(c)
                pos_diagonals.remove(r + c)
                neg_diagonals.remove(r - c)
        return False

    # Start placing queens from the first row
    # Ensure we skip rows that already have queens
    for r in range(n):
        if not any(q[0] == r for q in queens):
            if not place_queens(r):
                print(f"Failed to place a queen in row {r}")

    # Format the result
    result = ', '.join(f'{r} {c}' for r, c in sorted(queens))
    return f'<<<{result}>>>'

# Given board
board = [
    "000001X0",
    "00000000",
    "00000000",
    "X0000000",
    "00000000",
    "X0000000",
    "00000000",
    "000X0000"
]

# Solve and print the result
print(solve_n_queens(board))