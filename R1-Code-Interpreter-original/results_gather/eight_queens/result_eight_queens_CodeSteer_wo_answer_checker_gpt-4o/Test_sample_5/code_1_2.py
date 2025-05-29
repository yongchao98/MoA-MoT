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

def get_available_positions(board, row, n):
    return [col for col in range(n) if board[row][col] == 0 and is_safe(board, row, col, n)]

def solve_n_queens(board, n, result):
    if len(result) == n:
        return True

    # MRV: Choose the row with the fewest available positions
    row = min(range(n), key=lambda r: len(get_available_positions(board, r, n)) if r not in [r for r, _ in result] else float('inf'))

    # LCV: Try placing a queen in the least constraining column
    for col in sorted(get_available_positions(board, row, n), key=lambda c: sum(is_safe(board, r, c, n) for r in range(n))):
        board[row][col] = 1
        result.append((row, col))
        if solve_n_queens(board, n, result):
            return True
        board[row][col] = 0
        result.pop()

    return False

def place_queens():
    board = [
        [0, 0, 0, 0, 1, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, -1, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0]
    ]

    n = 8
    result = [(0, 4)]  # Initial queen position

    # Start solving
    if solve_n_queens(board, n, result):
        print("<<<" + ", ".join(f"{r} {c}" for r, c in result) + ">>>")
    else:
        print("No solution found")

place_queens()