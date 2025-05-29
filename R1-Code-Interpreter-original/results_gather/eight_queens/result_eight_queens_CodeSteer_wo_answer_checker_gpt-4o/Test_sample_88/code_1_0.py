def is_safe(board, row, col, n):
    # Check this column on upper side
    for i in range(row):
        if board[i][col] == 1:
            return False

    # Check upper diagonal on left side
    for i, j in zip(range(row, -1, -1), range(col, -1, -1)):
        if board[i][j] == 1:
            return False

    # Check upper diagonal on right side
    for i, j in zip(range(row, -1, -1), range(col, n)):
        if board[i][j] == 1:
            return False

    return True

def solve_n_queens(board, row, n, result):
    if row >= n:
        return True

    # Heuristic: Try columns with fewer constraints first
    columns = sorted(range(n), key=lambda c: sum(board[r][c] == 1 for r in range(n)))

    for col in columns:
        if board[row][col] == 0 and is_safe(board, row, col, n):
            board[row][col] = 1
            result.append((row, col))
            if solve_n_queens(board, row + 1, n, result):
                return True
            board[row][col] = 0
            result.pop()

    return False

def main():
    # Input board configuration with 'X' positions marked as -1
    board = [
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 1],
        [0, 0, 0, 0, 0, 0, 0, 0]
    ]

    # Example of marking 'X' positions (adjust as needed)
    # board[0][1] = -1  # Example of marking a position with 'X'
    # board[2][3] = -1  # Another example

    n = 8
    result = []

    if solve_n_queens(board, 0, n, result):
        print("<<<" + ", ".join(f"{r} {c}" for r, c in result) + ">>>")
    else:
        print("No solution found")

main()