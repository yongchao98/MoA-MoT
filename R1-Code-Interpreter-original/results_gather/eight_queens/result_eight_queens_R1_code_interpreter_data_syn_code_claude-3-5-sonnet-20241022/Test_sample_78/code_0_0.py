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
    # Initialize board with the given queen at (2,6)
    board = [[0 for x in range(n)] for y in range(n)]
    board[2][6] = 1  # Given queen
    
    def solve_util(queens_placed):
        if queens_placed == n:
            return True
        
        # If we already checked this position and placed a queen, skip it
        if queens_placed == 1 and board[2][6] == 1:
            return solve_util(queens_placed + 1)
        
        for row in range(n):
            for col in range(n):
                if board[row][col] == 0 and is_safe(board, row, col, n):
                    board[row][col] = 1
                    if solve_util(queens_placed + 1):
                        return True
                    board[row][col] = 0
        return False
    
    # Start with 1 queen already placed
    if solve_util(1):
        # Get positions of all queens
        positions = []
        for i in range(n):
            for j in range(n):
                if board[i][j] == 1:
                    positions.append(f"{i} {j}")
        print(",".join(positions))
    else:
        print("No solution exists")

solve_n_queens()