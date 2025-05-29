def solve_n_queens(board):
    n = len(board)
    queens = []
    cols = set()
    diag1 = set()
    diag2 = set()

    def is_valid(row, col):
        return col not in cols and (row - col) not in diag1 and (row + col) not in diag2 and board[row][col] != 'X'

    def place_queens(row):
        if row == n:
            return True
        for col in range(n):
            if is_valid(row, col):
                board[row][col] = 1
                queens.append((row, col))
                cols.add(col)
                diag1.add(row - col)
                diag2.add(row + col)
                if place_queens(row + 1):
                    return True
                # Backtrack
                board[row][col] = 0
                queens.pop()
                cols.remove(col)
                diag1.remove(row - col)
                diag2.remove(row + col)
        return False

    # Convert input to a mutable board
    board = [[c for c in line] for line in board]
    place_queens(0)

    # Collect the positions of the queens
    positions = [f"{r} {c}" for r, c in queens]
    return f"<<<{', '.join(positions)}>>>"

# Initial board setup
board = [
    "00100000",
    "00000000",
    "00000000",
    "0X000000",
    "000X0000",
    "00000000",
    "00000000",
    "000000X0"
]

# Solve the problem and print the result
print(solve_n_queens(board))