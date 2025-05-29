def is_safe(row, col, cols, pos_diags, neg_diags):
    return col not in cols and (row + col) not in pos_diags and (row - col) not in neg_diags

def solve_n_queens(board, row, positions, cols, pos_diags, neg_diags):
    if row >= len(board):
        return True

    for col in range(len(board)):
        if board[row][col] == 0 and is_safe(row, col, cols, pos_diags, neg_diags):
            board[row][col] = 1
            positions.append((row, col))
            cols.add(col)
            pos_diags.add(row + col)
            neg_diags.add(row - col)

            if solve_n_queens(board, row + 1, positions, cols, pos_diags, neg_diags):
                return True

            # Backtrack
            board[row][col] = 0
            positions.pop()
            cols.remove(col)
            pos_diags.remove(row + col)
            neg_diags.remove(row - col)

    return False

def place_queens(initial_board):
    board = [[0 if cell == '0' else cell for cell in row] for row in initial_board]
    positions = [(i, j) for i in range(8) for j in range(8) if board[i][j] == 1]

    # Initialize constraints
    cols = set(j for i, j in positions)
    pos_diags = set(i + j for i, j in positions)
    neg_diags = set(i - j for i, j in positions)

    if solve_n_queens(board, 0, positions, cols, pos_diags, neg_diags):
        result = ', '.join(f'{r} {c}' for r, c in positions)
        print(f'<<<{result}>>>')
    else:
        print("No solution exists")

# Example input with pre-existing queens and 'X' constraints
initial_board = [
    ['0', '0', '0', '0', '0', '0', '0', '0'],
    ['0', '0', '0', '0', '0', '0', '0', '0'],
    ['0', '0', '0', '0', '0', '0', '0', '0'],
    ['0', '0', '0', '0', '0', '0', '0', '0'],
    ['0', '0', '0', '0', '0', '0', '0', '0'],
    ['0', '0', '0', '0', '0', '0', '0', '0'],
    ['0', '0', '0', '0', '0', '0', '0', '0'],
    ['0', '0', '0', '0', '0', '0', '0', '0']
]

place_queens(initial_board)