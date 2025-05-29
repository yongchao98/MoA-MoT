def is_safe(board, row, col, n):
    # Check row
    for j in range(n):
        if board[row][j] == 1:
            return False
    
    # Check column
    for i in range(n):
        if board[i][col] == 1:
            return False
    
    # Check upper left diagonal
    for i, j in zip(range(row, -1, -1), range(col, -1, -1)):
        if board[i][j] == 1:
            return False
    
    # Check upper right diagonal
    for i, j in zip(range(row, -1, -1), range(col, n)):
        if board[i][j] == 1:
            return False
    
    # Check lower left diagonal
    for i, j in zip(range(row, n), range(col, -1, -1)):
        if board[i][j] == 1:
            return False
    
    # Check lower right diagonal
    for i, j in zip(range(row, n), range(col, n)):
        if board[i][j] == 1:
            return False
    
    return True

def solve_n_queens(board, n):
    # Count existing queens
    queen_count = sum(sum(row) for row in board)
    if queen_count == 8:
        return True
    
    # Try placing queens in remaining positions
    for row in range(n):
        for col in range(n):
            if board[row][col] == 0 and is_safe(board, row, col, n):
                board[row][col] = 1
                if solve_n_queens(board, n):
                    return True
                board[row][col] = 0
    
    return False

# Initialize the board with given queens
board = [
    [0,0,0,0,0,0,0,0],
    [0,0,1,0,0,0,0,0],
    [0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0],
    [0,0,0,0,1,0,0,0],
    [1,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0]
]

if solve_n_queens(board, 8):
    # Print queen positions
    positions = []
    for i in range(8):
        for j in range(8):
            if board[i][j] == 1:
                positions.append(f"{i} {j}")
    print(", ".join(positions))
else:
    print("No solution exists")