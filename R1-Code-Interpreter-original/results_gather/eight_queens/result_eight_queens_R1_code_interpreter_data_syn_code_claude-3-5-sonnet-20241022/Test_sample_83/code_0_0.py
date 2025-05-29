def is_safe(board, row, col, queens):
    # Check row
    for j in range(8):
        if board[row][j] == 1:
            return False
    
    # Check column
    for i in range(8):
        if board[i][col] == 1:
            return False
    
    # Check diagonal (top-left to bottom-right)
    for i, j in zip(range(row, -1, -1), range(col, -1, -1)):
        if board[i][j] == 1:
            return False
    for i, j in zip(range(row, 8), range(col, 8)):
        if board[i][j] == 1:
            return False
    
    # Check diagonal (top-right to bottom-left)
    for i, j in zip(range(row, -1, -1), range(col, 8)):
        if board[i][j] == 1:
            return False
    for i, j in zip(range(row, 8), range(col, -1, -1)):
        if board[i][j] == 1:
            return False
    
    return True

def solve_queens(board, queens_placed, queens):
    if queens_placed == 8:
        return True
    
    for row in range(8):
        for col in range(8):
            if board[row][col] == 0 and is_safe(board, row, col, queens):
                board[row][col] = 1
                queens.append((row, col))
                if solve_queens(board, queens_placed + 1, queens):
                    return True
                board[row][col] = 0
                queens.pop()
    
    return False

# Initialize the board with given queens
board = [
    [0,0,0,0,1,0,0,0],
    [0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,1],
    [0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0]
]

queens = [(0,4), (4,7)]
if solve_queens(board, 2, queens):
    # Sort queens by row for consistent output
    queens.sort()
    result = ", ".join(f"{row} {col}" for row, col in queens)
    print(f"<<<{result}>>>")
else:
    print("No solution exists")