def is_safe(board, row, col, n):
    # Check if a queen can be placed on board[row][col]
    
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
    # Initialize board with given configuration
    board = [[0 for x in range(n)] for y in range(n)]
    
    # Set initial queen and restricted positions
    board[1][4] = 1  # Initial queen
    restricted = [(4,7), (5,1), (7,0)]  # X positions
    
    queens = [(1,4)]  # Keep track of queen positions
    
    def solve_util(board, col):
        if len(queens) == 8:
            return True
            
        for row in range(n):
            if (row,col) not in restricted and board[row][col] == 0:
                if is_safe(board, row, col, n):
                    board[row][col] = 1
                    queens.append((row,col))
                    
                    if solve_util(board, (col + 1) % n):
                        return True
                        
                    board[row][col] = 0
                    queens.pop()
        return False
    
    # Start from column 0
    col = 0
    while len(queens) < 8:
        if not solve_util(board, col):
            col = (col + 1) % n
        
    # Format the output
    result = ", ".join(f"{row} {col}" for row, col in sorted(queens))
    print(f"<<<{result}>>>")

solve_n_queens()