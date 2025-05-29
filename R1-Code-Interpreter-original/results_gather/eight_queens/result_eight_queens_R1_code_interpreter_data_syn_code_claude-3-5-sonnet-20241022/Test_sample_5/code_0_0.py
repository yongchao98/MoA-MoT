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
    queens_positions = []
    
    # First, find existing queens
    for i in range(n):
        for j in range(n):
            if board[i][j] == 1:
                queens_positions.append((i, j))
    
    def solve(row):
        if len(queens_positions) == 8:
            return True
        
        if row >= n:
            return False
        
        for col in range(n):
            # Skip if position is blocked or already has a queen
            if board[row][col] == 'X' or board[row][col] == 1:
                continue
            
            if is_safe(board, row, col, n):
                board[row][col] = 1
                queens_positions.append((row, col))
                
                if solve(row + 1):
                    return True
                
                board[row][col] = 0
                queens_positions.pop()
        
        return False
    
    # Start solving from row 0
    solve(0)
    
    # Sort positions and format output
    queens_positions.sort()
    result = ", ".join(f"{row} {col}" for row, col in queens_positions)
    print(f"<<<{result}>>>")

# Initialize the board
board = [
    [0,0,0,0,0,1,0,0],
    [0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0],
    [0,'X',0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0]
]

solve_n_queens(board, 8)