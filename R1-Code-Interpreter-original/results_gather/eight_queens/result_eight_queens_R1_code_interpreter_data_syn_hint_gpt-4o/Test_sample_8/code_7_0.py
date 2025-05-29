def solve_n_queens_with_constraints(board):
    def is_safe(row, col):
        return not (cols[col] or diag1[row - col] or diag2[row + col] or board[row][col] == 'X')

    def place_queen(row):
        if row == 8:
            return True
        for col in range(8):
            if is_safe(row, col):
                queens.append((row, col))
                cols[col] = diag1[row - col] = diag2[row + col] = True
                if place_queen(row + 1):
                    return True
                queens.pop()
                cols[col] = diag1[row - col] = diag2[row + col] = False
        return False

    queens = []
    cols = [False] * 8
    diag1 = [False] * 15
    diag2 = [False] * 15

    place_queen(0)
    return queens

# Initial board setup with constraints
board = [
    ['0', '0', '0', '0', '0', '0', '0', '0'],
    ['0', '0', '0', '0', '0', '0', '0', '0'],
    ['0', '0', 'X', '0', '0', '0', '0', '0'],
    ['0', '0', '0', '0', '0', '0', '0', '0'],
    ['0', '0', '0', '0', '0', '0', '0', '0'],
    ['0', '0', '0', '0', '0', '0', '0', '0'],
    ['0', '0', '0', '0', '0', '0', '0', '0'],
    ['0', '0', '0', '0', '0', '0', '0', '0']
]

queen_positions = solve_n_queens_with_constraints(board)
print(queen_positions)