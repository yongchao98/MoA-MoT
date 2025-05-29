def is_safe(row, col, columns, diag1, diag2):
    if col in columns or (row - col) in diag1 or (row + col) in diag2:
        return False
    return True

def solve_n_queens(board, row, columns, diag1, diag2):
    if row == len(board):
        return True

    for col in range(len(board)):
        if board[row][col] == -1:
            continue
        if is_safe(row, col, columns, diag1, diag2):
            board[row][col] = 1
            columns.add(col)
            diag1.add(row - col)
            diag2.add(row + col)

            if solve_n_queens(board, row + 1, columns, diag1, diag2):
                return True

            board[row][col] = 0
            columns.remove(col)
            diag1.remove(row - col)
            diag2.remove(row + col)

    return False

def find_queen_positions():
    board = [
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, -1, 0, 0],
        [0, -1, 0, 0, 0, 0, -1, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, -1, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0]
    ]

    columns = set()
    diag1 = set()
    diag2 = set()

    if solve_n_queens(board, 0, columns, diag1, diag2):
        positions = []
        for i in range(len(board)):
            for j in range(len(board[i])):
                if board[i][j] == 1:
                    positions.append(f"{i} {j}")
        return positions
    else:
        return []

queen_positions = find_queen_positions()
print(queen_positions)