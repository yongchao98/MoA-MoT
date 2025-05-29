def is_valid(board, row, col):
    # Check column
    for i in range(row):
        if board[i] == col:
            return False
    # Check diagonals
    for i in range(row):
        if abs(board[i] - col) == abs(i - row):
            return False
    return True

def solve_n_queens(board, row, restricted):
    if row == 8:
        return True
    for col in range(8):
        if (row, col) in restricted:
            continue
        if is_valid(board, row, col):
            board[row] = col
            if solve_n_queens(board, row + 1, restricted):
                return True
            board[row] = -1
    return False

def main():
    # Initial board setup
    board = [-1] * 8
    board[4] = 0  # Queen already placed at (4,0)
    restricted = {(7, 5)}  # Restricted position

    # Solve the problem
    if solve_n_queens(board, 0, restricted):
        result = [f"{i} {board[i]}" for i in range(8)]
        print(f"<<<{', '.join(result)}>>>")
    else:
        print("No solution found")

main()