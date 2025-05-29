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
    for i, j in zip(range(row-1, -1, -1), range(col-1, -1, -1)):
        if board[i][j] == 1:
            return False
    
    # Check upper right diagonal
    for i, j in zip(range(row-1, -1, -1), range(col+1, n)):
        if board[i][j] == 1:
            return False
    
    # Check lower left diagonal
    for i, j in zip(range(row+1, n), range(col-1, -1, -1)):
        if board[i][j] == 1:
            return False
    
    # Check lower right diagonal
    for i, j in zip(range(row+1, n), range(col+1, n)):
        if board[i][j] == 1:
            return False
    
    return True

def solve_queens(board, queens_left, n):
    if queens_left == 0:
        return True
    
    for i in range(n):
        for j in range(n):
            # Skip if cell is already occupied or forbidden
            if board[i][j] != 0:
                continue
            
            if is_safe(board, i, j, n):
                board[i][j] = 1
                if solve_queens(board, queens_left - 1, n):
                    return True
                board[i][j] = 0
    
    return False

# Initialize board
n = 8
board = [
    [0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,1,0],
    [0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0],
    [1,0,0,-1,0,0,0,-1],  # -1 represents forbidden positions
    [0,0,0,0,0,0,0,0]
]

# Count existing queens
existing_queens = sum(row.count(1) for row in board)
queens_needed = 8 - existing_queens

if solve_queens(board, queens_needed, n):
    # Get positions of all queens
    positions = []
    for i in range(n):
        for j in range(n):
            if board[i][j] == 1:
                positions.append(f"{i} {j}")
    print(", ".join(positions))
else:
    print("No solution exists")