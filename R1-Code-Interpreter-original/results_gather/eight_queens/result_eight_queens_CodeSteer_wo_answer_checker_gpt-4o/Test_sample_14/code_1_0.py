def is_valid(board, row, col):
    # Check this column for any other queens
    for i in range(row):
        if board[i][col] == 1:
            return False

    # Check upper left diagonal for any other queens
    for i, j in zip(range(row, -1, -1), range(col, -1, -1)):
        if board[i][j] == 1:
            return False

    # Check upper right diagonal for any other queens
    for i, j in zip(range(row, -1, -1), range(col, len(board))):
        if board[i][j] == 1:
            return False

    return True

def solve_n_queens(board, row, positions):
    if row == len(board):
        return True

    for col in range(len(board)):
        # Check if the position is not marked with 'X' and is valid for placing a queen
        if board[row][col] == 0 and is_valid(board, row, col):
            board[row][col] = 1
            positions.append((row, col))
            if solve_n_queens(board, row + 1, positions):
                return True
            # Backtrack if placing a queen here doesn't lead to a solution
            board[row][col] = 0
            positions.pop()

    return False

def main():
    board = [
        [0, 0, 0, 0, 0, 0, 0, 0],
        [-1, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, -1, 0, 0, 0, 0],
        [0, -1, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0]
    ]

    positions = []
    if solve_n_queens(board, 0, positions):
        print("<<<" + ", ".join(f"{r} {c}" for r, c in positions) + ">>>")
    else:
        print("No solution found")

main()