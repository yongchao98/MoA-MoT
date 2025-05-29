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
    board = [[0 for x in range(n)] for y in range(n)]
    # Place the pre-existing queen
    board[7][5] = 1
    queens = [(7,5)]
    
    def solve_util(board, col, queens):
        if len(queens) == 8:
            return True
            
        for i in range(n):
            if is_safe(board, i, col, n):
                board[i][col] = 1
                queens.append((i,col))
                
                if solve_util(board, (col + 1) % n, queens):
                    return True
                    
                board[i][col] = 0
                queens.pop()
        
        return False
    
    # Start from column 0
    for start_col in range(n):
        if solve_util(board, start_col, queens):
            # Sort queens by row for consistent output
            queens.sort()
            result = ", ".join(f"{row} {col}" for row, col in queens)
            print(f"<<{result}>>")
            return
            
    print("No solution exists")

solve_n_queens()