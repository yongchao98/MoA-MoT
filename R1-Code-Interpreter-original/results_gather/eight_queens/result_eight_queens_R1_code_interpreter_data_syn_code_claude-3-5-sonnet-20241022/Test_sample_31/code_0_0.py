def is_safe(board, row, col, n):
    # Check row on left side
    for j in range(col):
        if board[row][j] == 1:
            return False
    
    # Check upper diagonal on left side
    for i, j in zip(range(row-1, -1, -1), range(col-1, -1, -1)):
        if board[i][j] == 1:
            return False
    
    # Check lower diagonal on left side
    for i, j in zip(range(row+1, n), range(col-1, -1, -1)):
        if board[i][j] == 1:
            return False
    
    return True

def solve_n_queens(n=8):
    # Initialize board with existing queens and constraints
    board = [[0 for x in range(n)] for y in range(n)]
    
    # Add constraints (X positions)
    board[0][5] = 2  # X position
    board[5][1] = 2  # X position
    
    # Add existing queen
    board[6][1] = 1
    
    def solve_util(board, col):
        if col >= n:
            return True
        
        # Skip column 1 as it already has a queen at row 6
        if col == 1:
            return solve_util(board, col + 1)
        
        for i in range(n):
            if board[i][col] != 2 and is_safe(board, i, col, n):  # Check if not X position
                board[i][col] = 1
                
                if solve_util(board, col + 1):
                    return True
                
                board[i][col] = 0
        
        return False
    
    if solve_util(board, 0):
        # Collect queen positions
        positions = []
        for i in range(n):
            for j in range(n):
                if board[i][j] == 1:
                    positions.append(f"{i} {j}")
        print(",".join(positions))
    else:
        print("Solution does not exist")

solve_n_queens()