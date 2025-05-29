def is_valid(board, row, col, n, cols, diag1, diag2):
    if cols[col] or diag1[row - col] or diag2[row + col]:
        return False
    return True

def solve_n_queens(board, row, n, positions, cols, diag1, diag2):
    if row == n:
        return len(positions) == n

    for col in range(n):
        if board[row][col] == 'X' or board[row][col] == 1:
            continue
        if is_valid(board, row, col, n, cols, diag1, diag2):
            board[row][col] = 1
            positions.append((row, col))
            cols[col] = diag1[row - col] = diag2[row + col] = True
            if solve_n_queens(board, row + 1, n, positions, cols, diag1, diag2):
                return True
            board[row][col] = 0
            positions.pop()
            cols[col] = diag1[row - col] = diag2[row + col] = False

    return False

def place_queens():
    board = [
        ['0', '0', 'X', '0', '0', '0', '0', '0'],
        ['1', '0', '0', '0', '0', '0', '0', '0'],
        ['0', '0', '0', '0', '0', '0', 'X', '0'],
        ['0', '0', '0', '0', '0', '0', '0', '0'],
        ['0', '0', '0', '0', '0', '0', '0', '0'],
        ['0', '0', '0', '0', '0', '0', '0', '0'],
        ['0', '0', '0', '0', '0', '0', '0', '0'],
        ['0', '0', '0', '0', '0', '0', '0', '0']
    ]

    n = 8
    positions = []

    # Initialize constraints
    cols = [False] * n
    diag1 = [False] * (2 * n)
    diag2 = [False] * (2 * n)

    # Add existing queens to positions and mark constraints
    for i in range(n):
        for j in range(n):
            if board[i][j] == '1':
                positions.append((i, j))
                cols[j] = True
                diag1[i - j] = True
                diag2[i + j] = True

    # Start solving from the first row
    solve_n_queens(board, 0, n, positions, cols, diag1, diag2)

    # Format the output
    result = ', '.join(f"{r} {c}" for r, c in positions)
    print(f"<<<{result}>>>")

place_queens()