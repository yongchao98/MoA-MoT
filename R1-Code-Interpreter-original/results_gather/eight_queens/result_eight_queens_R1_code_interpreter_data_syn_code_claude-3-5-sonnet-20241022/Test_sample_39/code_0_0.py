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
    # Initialize board with given constraints
    board = [
        [0,0,0,0,0,1,'X',0],
        [0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0],
        ['X',0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0],
        ['X',0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0],
        [0,0,0,'X',0,0,0,0]
    ]
    
    def solve_util(board, col):
        if col >= n:
            return True
            
        # Skip column 5 as queen is already placed there
        if col == 5:
            return solve_util(board, col + 1)
            
        for i in range(n):
            # Skip if position has 'X' or already has a queen
            if board[i][col] == 'X' or board[i][col] == 1:
                continue
                
            if is_safe(board, i, col, n):
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
        print("No solution exists")

solve_n_queens()