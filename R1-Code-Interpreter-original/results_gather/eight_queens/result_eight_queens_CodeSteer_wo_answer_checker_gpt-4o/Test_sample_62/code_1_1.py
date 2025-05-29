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

def solve_n_queens(board, row, n):
    if row >= n:
        return True

    for col in range(n):
        if board[row][col] == 0 and is_safe(board, row, col, n):
            board[row][col] = 1
            if solve_n_queens(board, row + 1, n):
                return True
            board[row][col] = 0

    return False

def place_queens():
    n = 8
    board = [[0 for _ in range(n)] for _ in range(n)]

    # Example of pre-placed queens and forbidden cells
    pre_placed_queens = [(0, 0), (1, 2)]  # Add tuples of (row, col) for pre-placed queens
    forbidden_cells = [(2, 3), (4, 5)]  # Add tuples of (row, col) for forbidden cells

    for row, col in pre_placed_queens:
        board[row][col] = 1

    for row, col in forbidden_cells:
        board[row][col] = 'X'

    if solve_n_queens(board, 0, n):
        result = []
        for i in range(n):
            for j in range(n):
                if board[i][j] == 1:
                    result.append(f"{i} {j}")
        print(f"<<<{', '.join(result)}>>>")
    else:
        print("No solution exists")

place_queens()