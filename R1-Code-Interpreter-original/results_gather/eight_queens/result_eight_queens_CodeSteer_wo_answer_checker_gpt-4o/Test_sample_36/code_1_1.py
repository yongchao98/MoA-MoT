def is_valid(board, row, col, queens):
    # Check if the column is safe
    for r in range(row):
        if queens[r] == col:
            return False

    # Check the diagonals
    for r in range(row):
        if abs(queens[r] - col) == abs(r - row):
            return False

    # Check if the position is marked with 'X'
    if board[row][col] == 'X':
        return False

    return True

def solve_n_queens(board, row, queens):
    if row == 8:
        return True

    for col in range(8):
        if is_valid(board, row, col, queens):
            queens[row] = col
            if solve_n_queens(board, row + 1, queens):
                return True
            queens[row] = -1

    return False

def main():
    board = [
        "00000000",
        "00000000",
        "000X0000",
        "0000000X",
        "00000000",
        "0000000X",
        "00000000",
        "0000X000"
    ]

    queens = [-1] * 8
    if solve_n_queens(board, 0, queens):
        result = ', '.join(f"{r} {c}" for r, c in enumerate(queens))
        print(f"<<<{result}>>>")
    else:
        print("No solution found")

main()