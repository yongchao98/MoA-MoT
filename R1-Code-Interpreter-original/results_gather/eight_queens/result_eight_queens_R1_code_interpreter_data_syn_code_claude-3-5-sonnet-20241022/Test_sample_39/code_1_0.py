def is_safe(board, row, col, n):
    # Check row
    for j in range(n):
        if board[row][j] == 1:
            return False
    
    # Check column
    for i in range(n):
        if board[i][col] == 1:
            return False
    
    # Check upper right diagonal
    for i, j in zip(range(row-1, -1, -1), range(col+1, n)):
        if i >= 0 and j < n and board[i][j] == 1:
            return False
    
    # Check upper left diagonal
    for i, j in zip(range(row-1, -1, -1), range(col-1, -1, -1)):
        if i >= 0 and j >= 0 and board[i][j] == 1:
            return False
    
    # Check lower right diagonal
    for i, j in zip(range(row+1, n), range(col+1, n)):
        if i < n and j < n and board[i][j] == 1:
            return False
    
    # Check lower left diagonal
    for i, j in zip(range(row+1, n), range(col-1, -1, -1)):
        if i < n and j >= 0 and board[i][j] == 1:
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
    
    queens_placed = 1  # Start with 1 as we already have one queen
    
    def solve_util(board, queens_placed):
        if queens_placed == 8:
            return True
            
        for i in range(n):
            for j in range(n):
                # Skip if position has 'X' or is the pre-placed queen position
                if board[i][j] == 'X' or (i == 0 and j == 5):
                    continue
                    
                if board[i][j] == 0 and is_safe(board, i, j, n):
                    board[i][j] = 1
                    
                    if solve_util(board, queens_placed + 1):
                        return True
                        
                    board[i][j] = 0
        
        return False
    
    if solve_util(board, queens_placed):
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