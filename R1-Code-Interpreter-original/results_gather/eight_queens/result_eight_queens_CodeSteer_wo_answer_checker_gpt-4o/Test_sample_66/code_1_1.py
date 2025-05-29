def is_valid(board, row, col, queens):
    # Check if the column is safe
    for r, c in queens:
        if c == col or abs(r - row) == abs(c - col):
            return False
    return True

def solve_n_queens(board, row, queens):
    if row == 8:
        return True
    
    for col in range(8):
        if board[row][col] == 0 and is_valid(board, row, col, queens):
            queens.append((row, col))
            if solve_n_queens(board, row + 1, queens):
                return True
            queens.pop()
    
    return False

def main():
    # Initialize the board with the 'X' position
    board = [
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 'X', 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0]
    ]
    
    queens = []
    if solve_n_queens(board, 0, queens):
        result = ', '.join(f"{r} {c}" for r, c in queens)
        print(f"<<<{result}>>>")
    else:
        print("No solution found")

main()