def solve_n_queens(board):
    def is_safe(row, col):
        return col not in cols and (row - col) not in neg_diags and (row + col) not in pos_diags and board[row][col] != 'X'

    def place_queens(row):
        if row == 8:
            return True
        for col in range(8):
            if is_safe(row, col):
                queens.append((row, col))
                cols.add(col)
                neg_diags.add(row - col)
                pos_diags.add(row + col)
                if place_queens(row + 1):
                    return True
                queens.pop()
                cols.remove(col)
                neg_diags.remove(row - col)
                pos_diags.remove(row + col)
        return False

    queens = []
    cols = set()
    neg_diags = set()
    pos_diags = set()

    # Mark the initial board constraints
    for r in range(8):
        for c in range(8):
            if board[r][c] == '1':
                queens.append((r, c))
                cols.add(c)
                neg_diags.add(r - c)
                pos_diags.add(r + c)

    # Start placing queens from the first row
    place_queens(0)

    # Format the output
    return '<<<' + ', '.join(f'{r} {c}' for r, c in queens) + '>>>'

# Initial board setup
board = [
    ['0', '0', '0', '0', '0', '0', '0', '0'],
    ['0', '0', '0', '0', '0', '0', '0', '0'],
    ['0', '0', '0', '0', '0', '0', '0', '0'],
    ['0', '0', '0', '0', '0', '0', '0', '0'],
    ['0', '0', '0', '0', '0', '0', '0', '0'],
    ['X', 'X', 'X', '0', '0', '0', '0', '0'],
    ['0', '0', '0', '0', '0', '0', '0', '0'],
    ['0', '0', '0', '0', '0', '0', '0', '0']
]

# Solve the problem
print(solve_n_queens(board))