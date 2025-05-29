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
    
    def solve_util(board, col):
        if len(queens) == 8:
            return True
            
        for i in range(n):
            # Skip row 7 as it already has a queen
            if i == 7:
                continue
            if is_safe(board, i, col, n):
                board[i][col] = 1
                queens.append((i,col))
                
                next_col = col + 1
                if next_col == 5:  # Skip column 5 as it already has a queen
                    next_col += 1
                if next_col < n and solve_util(board, next_col):
                    return True
                    
                board[i][col] = 0
                queens.pop()
        
        return False
    
    # Start from column 0
    if solve_util(board, 0):
        # Sort queens by row for consistent output
        queens.sort()
        result = ", ".join(f"{row} {col}" for row, col in queens)
        print(f"<<{result}>>")
    else:
        print("No solution exists")

solve_n_queens()