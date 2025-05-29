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
    
    # Check column above
    for i in range(row):
        if board[i][col] == 1:
            return False
            
    # Check column below
    for i in range(row+1, n):
        if board[i][col] == 1:
            return False
    
    # Check upper diagonal on right side
    for i, j in zip(range(row-1, -1, -1), range(col+1, n)):
        if board[i][j] == 1:
            return False
    
    # Check lower diagonal on right side
    for i, j in zip(range(row+1, n), range(col+1, n)):
        if board[i][j] == 1:
            return False
    
    return True

def solve_n_queens(n=8):
    # Initialize board with given queen
    board = [[0 for x in range(n)] for y in range(n)]
    board[0][4] = 1  # Initial queen position
    
    def solve_util(board, col):
        if col >= n:
            return True
            
        # Skip column 4 in first row as it already has a queen
        if col == 4:
            return solve_util(board, col + 1)
            
        for i in range(n):
            # Skip first row as it already has a queen
            if col < n and i == 0:
                continue
                
            if is_safe(board, i, col, n):
                board[i][col] = 1
                
                if solve_util(board, col + 1):
                    return True
                    
                board[i][col] = 0
        
        return False
    
    if solve_util(board, 0):
        # Convert board to position list
        positions = []
        for i in range(n):
            for j in range(n):
                if board[i][j] == 1:
                    positions.append(f"{i} {j}")
        print(",".join(positions))
    else:
        print("No solution exists")

solve_n_queens()