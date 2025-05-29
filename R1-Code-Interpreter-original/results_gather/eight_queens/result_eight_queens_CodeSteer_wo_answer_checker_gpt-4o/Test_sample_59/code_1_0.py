def is_safe(row, col, cols, diag1, diag2):
    return col not in cols and (row - col) not in diag1 and (row + col) not in diag2

def solve_n_queens(board, row, cols, diag1, diag2, positions):
    if row == len(board):
        return True

    for col in range(len(board)):
        if board[row][col] == 0 and is_safe(row, col, cols, diag1, diag2):
            board[row][col] = 1
            positions.append((row, col))
            cols.add(col)
            diag1.add(row - col)
            diag2.add(row + col)

            if solve_n_queens(board, row + 1, cols, diag1, diag2, positions):
                return True

            # Backtrack
            board[row][col] = 0
            positions.pop()
            cols.remove(col)
            diag1.remove(row - col)
            diag2.remove(row + col)

    return False

def place_queens_with_constraints(initial_board):
    board = [[0 if cell == '0' else -1 for cell in row] for row in initial_board]
    positions = []
    cols = set()
    diag1 = set()
    diag2 = set()

    # Place initial queens and mark their attacks
    for r in range(8):
        for c in range(8):
            if board[r][c] == 1:
                positions.append((r, c))
                cols.add(c)
                diag1.add(r - c)
                diag2.add(r + c)

    if solve_n_queens(board, 0, cols, diag1, diag2, positions):
        result = ', '.join(f'{r} {c}' for r, c in positions)
        print(f'<<<{result}>>>')
    else:
        print("No solution exists")

# Example initial board with no queens or 'X' positions
initial_board = [
    "00000000",
    "00000000",
    "00000000",
    "00000000",
    "00000000",
    "00000000",
    "00000000",
    "00000000"
]

place_queens_with_constraints(initial_board)