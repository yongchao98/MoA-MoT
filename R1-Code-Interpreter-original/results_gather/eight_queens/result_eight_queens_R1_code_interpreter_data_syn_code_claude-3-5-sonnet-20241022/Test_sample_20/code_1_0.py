def is_safe(board, row, col, n):
    # Check row
    for j in range(n):
        if board[row][j] == 1:
            return False
    
    # Check column
    for i in range(n):
        if board[i][col] == 1:
            return False
    
    # Check upper diagonal
    for i, j in zip(range(row, -1, -1), range(col, -1, -1)):
        if board[i][j] == 1:
            return False
    
    # Check lower diagonal
    for i, j in zip(range(row, n), range(col, -1, -1)):
        if board[i][j] == 1:
            return False
            
    # Check upper diagonal (right side)
    for i, j in zip(range(row, -1, -1), range(col, n)):
        if board[i][j] == 1:
            return False
            
    # Check lower diagonal (right side)
    for i, j in zip(range(row, n), range(col, n)):
        if board[i][j] == 1:
            return False
    
    return True

def solve_n_queens(n=8):
    board = [[0]*n for _ in range(n)]
    # Place initial queen
    board[0][1] = 1
    # Mark forbidden position
    forbidden_row, forbidden_col = 1, 5
    
    queens_placed = 1
    
    def solve_util(board, queens_placed):
        if queens_placed == 8:
            return True
            
        for i in range(n):
            for j in range(n):
                # Skip if position is initial queen position or forbidden position
                if (i == 0 and j == 1) or (i == forbidden_row and j == forbidden_col):
                    continue
                    
                if board[i][j] == 0 and is_safe(board, i, j, n):
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