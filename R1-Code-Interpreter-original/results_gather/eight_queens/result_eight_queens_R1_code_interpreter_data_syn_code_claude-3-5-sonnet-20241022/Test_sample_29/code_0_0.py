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

def solve_n_queens(board, n):
    queens_positions = []
    # First, find existing queens
    for i in range(n):
        for j in range(n):
            if board[i][j] == 1:
                queens_positions.append((i, j))
    
    def solve(row, queens_left):
        if queens_left == 0:
            return True
        
        if row >= n:
            return False
            
        # If this row already has a queen, skip to next row
        for pos in queens_positions:
            if pos[0] == row:
                return solve(row + 1, queens_left)
        
        for col in range(n):
            if is_safe(board, row, col, n):
                board[row][col] = 1
                queens_positions.append((row, col))
                if solve(row + 1, queens_left - 1):
                    return True
                board[row][col] = 0
                queens_positions.pop()
        
        return solve(row + 1, queens_left)
    
    # Start solving from row 0 with 6 queens left to place
    solve(0, 6)
    
    # Sort positions and format output
    queens_positions.sort()
    result = ", ".join(f"{pos[0]} {pos[1]}" for pos in queens_positions)
    print(f"<<<{result}>>>")

# Initialize the board
board = [
    [0,0,0,0,0,0,0,0],
    [1,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0],
    [0,0,0,0,0,1,0,0],
    [0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0],
    [0,1,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0]
]

solve_n_queens(board, 8)