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
    
    # Check row on right side
    for j in range(col+1, n):
        if board[row][j] == 1:
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
    # Initialize board with pre-placed queen and blocked cell
    board = [[0 for x in range(n)] for y in range(n)]
    board[4][2] = 1  # Pre-placed queen
    blocked = (5, 5)  # Blocked cell
    
    def solve_util(col):
        if col >= n:
            return True
        
        # Skip column 2 as it already has a queen
        if col == 2:
            return solve_util(col + 1)
        
        for i in range(n):
            # Skip if position is blocked
            if (i, col) == blocked:
                continue
                
            if is_safe(board, i, col, n):
                board[i][col] = 1
                
                if solve_util(col + 1):
                    return True
                
                board[i][col] = 0
        
        return False
    
    if solve_util(0):
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