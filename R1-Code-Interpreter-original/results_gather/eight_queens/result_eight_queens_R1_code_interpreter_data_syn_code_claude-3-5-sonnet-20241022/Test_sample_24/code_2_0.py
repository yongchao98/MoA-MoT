def is_safe(board, row, col):
    n = 8
    # Check if a queen can be placed on board[row][col]
    
    # Check this row on left side
    for j in range(col):
        if board[row][j] == 1:
            return False
            
    # Check this row on right side
    for j in range(col + 1, n):
        if board[row][j] == 1:
            return False
    
    # Check this column
    for i in range(n):
        if i != row and board[i][col] == 1:
            return False
            
    # Check upper diagonal on left side
    i, j = row - 1, col - 1
    while i >= 0 and j >= 0:
        if board[i][j] == 1:
            return False
        i -= 1
        j -= 1
    
    # Check lower diagonal on left side
    i, j = row + 1, col - 1
    while i < n and j >= 0:
        if board[i][j] == 1:
            return False
        i += 1
        j -= 1
        
    # Check upper diagonal on right side
    i, j = row - 1, col + 1
    while i >= 0 and j < n:
        if board[i][j] == 1:
            return False
        i -= 1
        j += 1
        
    # Check lower diagonal on right side
    i, j = row + 1, col + 1
    while i < n and j < n:
        if board[i][j] == 1:
            return False
        i += 1
        j += 1
    
    return True

def solve_queens():
    n = 8
    # Initialize board with given constraints
    board = [
        [0,0,0,-1,0,0,0,0],
        [0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0],
        [0,0,1,0,-1,-1,0,0],
        [0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0]
    ]
    
    queens_needed = 7  # We already have 1 queen
    
    def solve(queens_left):
        if queens_left == 0:
            return True
            
        # Try each cell on the board
        for i in range(n):
            for j in range(n):
                # Skip if cell is not empty or forbidden
                if board[i][j] != 0:
                    continue
                    
                if is_safe(board, i, j):
                    # Place queen
                    board[i][j] = 1
                    
                    if solve(queens_left - 1):
                        return True
                        
                    # If placing queen didn't lead to a solution, remove it
                    board[i][j] = 0
                    
        return False
    
    if solve(queens_needed):
        # Collect queen positions
        positions = []
        for i in range(n):
            for j in range(n):
                if board[i][j] == 1:
                    positions.append(f"{i} {j}")
        print(",".join(positions))
    else:
        print("No solution exists")

solve_queens()