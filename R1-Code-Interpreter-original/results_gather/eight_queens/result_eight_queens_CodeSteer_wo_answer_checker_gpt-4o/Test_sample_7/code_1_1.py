def is_valid(board, row, col):
    for i in range(row):
        if board[i] == col or abs(board[i] - col) == abs(i - row):
            return False
    return True

def get_mrv_row(board, constraints):
    min_options = float('inf')
    mrv_row = -1
    for row in range(8):
        if board[row] == -1:
            options = sum(constraints[row])
            if options < min_options:
                min_options = options
                mrv_row = row
    return mrv_row

def update_constraints(constraints, row, col, add=True):
    for r in range(8):
        if r != row:
            constraints[r][col] += 1 if add else -1
            if col - (row - r) >= 0:
                constraints[r][col - (row - r)] += 1 if add else -1
            if col + (row - r) < 8:
                constraints[r][col + (row - r)] += 1 if add else -1

def solve_n_queens(board, constraints, positions):
    if len(positions) == 8:
        return True
    row = get_mrv_row(board, constraints)
    if row == -1:
        return False
    for col in sorted(range(8), key=lambda c: constraints[row][c]):
        if (row, col) == (2, 5) or constraints[row][col] > 0:
            continue
        if is_valid(board, row, col):
            board[row] = col
            positions.append((row, col))
            update_constraints(constraints, row, col, add=True)
            if solve_n_queens(board, constraints, positions):
                return True
            update_constraints(constraints, row, col, add=False)
            positions.pop()
            board[row] = -1
    return False

def place_queens():
    board = [-1] * 8
    constraints = [[0] * 8 for _ in range(8)]
    positions = []
    solve_n_queens(board, constraints, positions)
    return positions

queen_positions = place_queens()
formatted_positions = ', '.join(f'{r} {c}' for r, c in queen_positions)
print(f'<<<{formatted_positions}>>>')