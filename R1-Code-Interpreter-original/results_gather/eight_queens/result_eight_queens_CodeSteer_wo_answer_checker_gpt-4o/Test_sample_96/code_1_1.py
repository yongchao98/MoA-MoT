def solve_n_queens(n, row, cols, diags1, diags2, positions):
    if row == n:
        return True

    for col in range(n):
        if (cols & (1 << col)) or (diags1 & (1 << (row - col + n - 1))) or (diags2 & (1 << (row + col))):
            continue

        # Place the queen
        positions.append((row, col))
        cols |= (1 << col)
        diags1 |= (1 << (row - col + n - 1))
        diags2 |= (1 << (row + col))

        if solve_n_queens(n, row + 1, cols, diags1, diags2, positions):
            return True

        # Remove the queen (backtrack)
        positions.pop()
        cols &= ~(1 << col)
        diags1 &= ~(1 << (row - col + n - 1))
        diags2 &= ~(1 << (row + col))

    return False

def place_queens():
    n = 8
    board = [
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 'X', 0, 0],
        [0, 'X', 0, 0, 0, 0, 'X', 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 'X', 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0]
    ]

    # Convert 'X' to -1 for easier processing
    for i in range(n):
        for j in range(n):
            if board[i][j] == 'X':
                board[i][j] = -1

    positions = []
    cols = diags1 = diags2 = 0

    if solve_n_queens(n, 0, cols, diags1, diags2, positions):
        result = ', '.join(f"{r} {c}" for r, c in positions)
        print(f"<<<{result}>>>")
    else:
        print("No solution found")

place_queens()