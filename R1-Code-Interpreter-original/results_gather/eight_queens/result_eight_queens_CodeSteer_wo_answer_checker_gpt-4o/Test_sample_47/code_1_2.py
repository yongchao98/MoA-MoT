def is_safe(board, row, col):
    # Check if the position is marked with 'X'
    if board[row][col] == 'X':
        return False

    # Check this column on upper side
    for i in range(row):
        if board[i][col] == 1:
            return False

    # Check upper diagonal on left side
    for i, j in zip(range(row, -1, -1), range(col, -1, -1)):
        if board[i][j] == 1:
            return False

    # Check upper diagonal on right side
    for i, j in zip(range(row, -1, -1), range(col, len(board), 1)):
        if board[i][j] == 1:
            return False

    return True

def solve_n_queens(board, row, queens_positions):
    if row >= len(board):
        return True

    for col in range(len(board)):
        if is_safe(board, row, col):
            board[row][col] = 1
            queens_positions.append((row, col))

            if solve_n_queens(board, row + 1, queens_positions):
                return True

            board[row][col] = 0
            queens_positions.pop()

    return False

def main():
    board = [
        [0, 0, 0, 0, 0, 0, 0, 'X'],
        [0, 0, 0, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        ['X', 0, 0, 'X', 0, 0, 0, 0],
        [0, 0, 0, 'X', 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0]
    ]

    queens_positions = [(1, 3)]
    board[1][3] = 1

    if solve_n_queens(board, 0, queens_positions):
        result = ', '.join(f"{r} {c}" for r, c in queens_positions)
        print(f"<<<{result}>>>")
    else:
        print("No solution found")

main()