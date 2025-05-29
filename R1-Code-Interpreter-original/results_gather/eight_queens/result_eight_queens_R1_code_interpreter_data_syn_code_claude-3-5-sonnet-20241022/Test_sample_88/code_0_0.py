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

def solve_queens(board, n):
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
            
        # If this row already has a queen, move to next row
        if any(pos[0] == row for pos in queens_positions):
            return solve(row + 1)
        
        for col in range(n):
            if is_safe(board, row, col, n):
                board[row][col] = 1
                queens_positions.append((row, col))
                if solve(row + 1):
                    return True
                board[row][col] = 0
                queens_positions.pop()
        return False
    
    solve(0)
    return queens_positions

# Initialize the board
board = [
    [0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0],
    [0,0,1,0,0,0,0,0],
    [0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,1],
    [0,0,0,0,0,0,0,0]
]

positions = solve_queens(board, 8)
# Sort positions by row for consistent output
positions.sort()
# Convert to string format
result = ", ".join(f"{pos[0]} {pos[1]}" for pos in positions)
print(result)