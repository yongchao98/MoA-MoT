def is_safe(board, row, col, n):
    # Check row on left and right side
    for j in range(n):
        if board[row][j] == 1:
            return False
    
    # Check column
    for i in range(n):
        if board[i][col] == 1:
            return False
    
    # Check diagonals
    for i in range(n):
        for j in range(n):
            if board[i][j] == 1:
                if abs(row - i) == abs(col - j):
                    return False
    
    return True

def solve_n_queens(n=8):
    # Initialize board with pre-placed queen and forbidden positions
    board = [
        [0,0,0,0,0,0,0,-1],  # X at (0,7)
        [0,0,0,1,0,0,0,0],   # Pre-placed queen at (1,3)
        [-1,0,0,0,0,0,0,0],  # X at (2,0)
        [0,0,0,0,0,0,0,0],
        [0,0,-1,0,0,0,0,0],  # X at (4,2)
        [0,0,0,0,-1,0,0,0],  # X at (5,4)
        [0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0]
    ]
    
    queens_placed = 1  # Start with 1 as we already have one queen
    
    def solve_util(board, queens_placed):
        if queens_placed == 8:
            return True
            
        for i in range(n):
            for j in range(n):
                # Skip if position is forbidden or already has a queen
                if board[i][j] == -1 or board[i][j] == 1:
                    continue
                    
                if is_safe(board, i, j, n):
                    board[i][j] = 1
                    
                    if solve_util(board, queens_placed + 1):
                        return True
                        
                    board[i][j] = 0
        
        return False
    
    if solve_util(board, queens_placed):
        positions = []
        for i in range(n):
            for j in range(n):
                if board[i][j] == 1:
                    positions.append(f"{i} {j}")
        print(",".join(positions))
    else:
        print("No solution exists")

solve_n_queens()