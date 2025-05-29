def solve_n_queens(board):
    def is_safe(row, col):
        # Check if the column is safe
        if col in cols:
            return False
        # Check if the main diagonal is safe
        if (row - col) in main_diag:
            return False
        # Check if the anti-diagonal is safe
        if (row + col) in anti_diag:
            return False
        return True

    def place_queens(row):
        if row == 8:
            return True
        for col in range(8):
            if board[row][col] == 'X' or board[row][col] == '1' or not is_safe(row, col):
                continue
            # Place the queen
            board[row][col] = '1'
            cols.add(col)
            main_diag.add(row - col)
            anti_diag.add(row + col)
            print(f"Placing queen at ({row}, {col})")
            # Recur to place the rest
            if place_queens(row + 1):
                return True
            # Backtrack
            print(f"Backtracking from ({row}, {col})")
            board[row][col] = '0'
            cols.remove(col)
            main_diag.remove(row - col)
            anti_diag.remove(row + col)
        return False

    # Initialize the board
    board = [list(row) for row in board]
    cols = set()
    main_diag = set()
    anti_diag = set()

    # Mark the existing queens
    for r in range(8):
        for c in range(8):
            if board[r][c] == '1':
                cols.add(c)
                main_diag.add(r - c)
                anti_diag.add(r + c)

    # Solve the problem
    if not place_queens(0):
        print("No solution found")

    # Collect the positions of the queens
    result = []
    for r in range(8):
        for c in range(8):
            if board[r][c] == '1':
                result.append(f"{r} {c}")

    return f"<<<{', '.join(result)}>>>"

# Initial board setup
initial_board = [
    "00000000",
    "00001000",
    "00000000",
    "X0000000",
    "00000000",
    "00000000",
    "00000000",
    "00000000"
]

# Solve and print the result
print(solve_n_queens(initial_board))