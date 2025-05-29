def is_safe(board, row, col, n):
    # Check row on left side
    for j in range(col):
        if board[row][j] == 1:
            return False
    
    # Check row on right side
    for j in range(col+1, n):
        if board[row][j] == 1:
            return False
            
    # Check column
    for i in range(n):
        if board[i][col] == 1:
            return False
    
    # Check upper diagonal on left side
    i, j = row-1, col-1
    while i >= 0 and j >= 0:
        if board[i][j] == 1:
            return False
        i -= 1
        j -= 1
    
    # Check lower diagonal on left side
    i, j = row+1, col-1
    while i < n and j >= 0:
        if board[i][j] == 1:
            return False
        i += 1
        j -= 1
    
    # Check upper diagonal on right side
    i, j = row-1, col+1
    while i >= 0 and j < n:
        if board[i][j] == 1:
            return False
        i -= 1
        j += 1
    
    # Check lower diagonal on right side
    i, j = row+1, col+1
    while i < n and j < n:
        if board[i][j] == 1:
            return False
        i += 1
        j += 1
    
    return True

def solve_n_queens(n=8):
    board = [[0 for x in range(n)] for y in range(n)]
    board[2][6] = 1  # Given queen
    queens_positions = [(2,6)]
    
    def solve_util(col):
        if len(queens_positions) == n:
            return True
            
        # Skip column 6 as it already has a queen at row 2
        if col == 6:
            return solve_util(col + 1)
            
        if col >= n:
            return False
            
        for row in range(n):
            if is_safe(board, row, col, n):
                board[row][col] = 1
                queens_positions.append((row,col))
                
                if solve_util(col + 1):
                    return True
                    
                board[row][col] = 0
                queens_positions.pop()
                
        return False
    
    if solve_util(0):
        # Sort positions by row for consistent output
        queens_positions.sort()
        print(",".join(f"{row} {col}" for row, col in queens_positions))
    else:
        print("No solution exists")

solve_n_queens()