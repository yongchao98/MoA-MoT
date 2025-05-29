def is_safe(board, row, col, n):
    # Check this row on left side
    for i in range(col):
        if board[row][i] == 1:
            return False

    # Check upper diagonal on left side
    for i, j in zip(range(row, -1, -1), range(col, -1, -1)):
        if board[i][j] == 1:
            return False

    # Check lower diagonal on left side
    for i, j in zip(range(row, n, 1), range(col, -1, -1)):
        if board[i][j] == 1:
            return False

    return True

def solve_n_queens_util(board, col, n, blocked_positions):
    # base case: If all queens are placed
    if col >= n:
        return True

    # If there's already a queen in this column, move to the next column
    for i in range(n):
        if board[i][col] == 1:
            return solve_n_queens_util(board, col + 1, n, blocked_positions)

    for i in range(n):
        if (i, col) not in blocked_positions and is_safe(board, i, col, n):
            # Place this queen in board[i][col]
            board[i][col] = 1

            # Recur to place rest of the queens
            if solve_n_queens_util(board, col + 1, n, blocked_positions):
                return True

            # If placing queen in board[i][col] doesn't lead to a solution
            # then remove queen from board[i][col]
            board[i][col] = 0

    # If the queen cannot be placed in any row in this column col then return false
    return False

def solve_n_queens(n, initial_board):
    board = [[0] * n for _ in range(n)]
    blocked_positions = set()
    initial_queens = set()
    for i in range(n):
        for j in range(n):
            if initial_board[i][j] == '1':
                board[i][j] = 1
                initial_queens.add((i, j))
            elif initial_board[i][j] == 'X':
                blocked_positions.add((i, j))

    # Start solving from the first column that doesn't have a queen
    start_col = 0
    while start_col < n and any(board[i][start_col] == 1 for i in range(n)):
        start_col += 1

    if not solve_n_queens_util(board, start_col, n, blocked_positions):
        return "No solution exists"

    # Collect the positions of the queens
    queen_positions = []
    for i in range(n):
        for j in range(n):
            if board[i][j] == 1:
                queen_positions.append(f"{i} {j}")

    return queen_positions

# Initial board setup
initial_board = [
    "000X0000",
    "00000000",
    "00000000",
    "00000000",
    "0010XX00",
    "00000000",
    "00000000",
    "00000000"
]

# Solve the problem
queen_positions = solve_n_queens(8, initial_board)
print(queen_positions)