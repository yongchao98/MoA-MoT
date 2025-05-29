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
        if board[i][j] == 1:
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
    # Initialize board with pre-placed queen and X
    board = [[0 for x in range(n)] for y in range(n)]
    board[1][4] = 1  # Pre-placed queen
    board[3][0] = 'X'  # Forbidden position
    
    def solve_util(board, col):
        if col >= n:
            return True
        
        # Skip column 4 as it already has a queen in row 1
        if col == 4:
            return solve_util(board, col + 1)
            
        for i in range(n):
            # Skip row 1 as it already has a queen
            if col != 4 and i == 1:
                continue
            # Skip forbidden position
            if i == 3 and col == 0:
                continue
                
            if is_safe(board, i, col, n) and board[i][col] != 'X':
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
        print(", ".join(positions))
    else:
        print("No solution exists")

solve_n_queens()