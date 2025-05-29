def is_valid(board, row, col, queens):
    # Check if there's a queen in the same column or diagonal
    for r, c in queens:
        if c == col or abs(r - row) == abs(c - col):
            return False
    return True

def solve_n_queens(board, row, queens):
    if row == 8:
        return queens
    
    # Skip rows that already have a queen
    if any(board[row][col] == '1' for col in range(8)):
        return solve_n_queens(board, row + 1, queens)
    
    for col in range(8):
        if board[row][col] == 'X' or board[row][col] == '1':
            continue
        if is_valid(board, row, col, queens):
            queens.append((row, col))
            result = solve_n_queens(board, row + 1, queens)
            if result:
                return result
            queens.pop()
    return None

def main():
    board = [
        "00000000",
        "00000000",
        "0X0000X0",
        "00000X00",
        "00000000",
        "00000000",
        "0000000X",
        "00000000"
    ]
    
    # Parse the board to find pre-placed queens
    queens = [(r, c) for r in range(8) for c in range(8) if board[r][c] == '1']
    
    # Start solving from the first row
    solution = solve_n_queens(board, 0, queens)
    if solution:
        print("<<<" + ", ".join(f"{r} {c}" for r, c in solution) + ">>>")
    else:
        print("No solution found")

main()