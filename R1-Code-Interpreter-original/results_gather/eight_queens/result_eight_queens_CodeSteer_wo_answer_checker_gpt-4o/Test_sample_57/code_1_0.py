def is_safe(row, col, cols, diag1, diag2, board):
    # Check if the column, and both diagonals are free and not marked 'X'
    return col not in cols and (row - col) not in diag1 and (row + col) not in diag2 and board[row][col] != 'X'

def solve_n_queens(board, row, positions, cols, diag1, diag2):
    if row == len(board):
        return True

    for col in range(len(board)):
        if is_safe(row, col, cols, diag1, diag2, board):
            # Place the queen
            board[row][col] = 1
            positions.append((row, col))
            cols.add(col)
            diag1.add(row - col)
            diag2.add(row + col)

            if solve_n_queens(board, row + 1, positions, cols, diag1, diag2):
                return True

            # Backtrack
            board[row][col] = 0
            positions.pop()
            cols.remove(col)
            diag1.remove(row - col)
            diag2.remove(row + col)

    return False

def place_queens():
    # Initialize the board with 0s and mark the 'X' position
    board = [[0] * 8 for _ in range(8)]
    board[5][3] = 'X'  # Mark the 'X' position

    # Initialize sets to track occupied columns and diagonals
    cols = set()
    diag1 = set()
    diag2 = set()

    # List to store the positions of the queens
    positions = []

    # Start solving from the first row
    solve_n_queens(board, 0, positions, cols, diag1, diag2)

    # Format the output
    result = ', '.join(f"{r} {c}" for r, c in positions)
    print(f"<<<{result}>>>")

place_queens()