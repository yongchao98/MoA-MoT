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
    
    return True

def solve_n_queens(n=8):
    board = [[0 for x in range(n)] for y in range(n)]
    # Mark the X position
    forbidden = (2, 2)
    
    def solve_util(board, col):
        if col >= n:
            return True
        
        for i in range(n):
            if (i, col) == forbidden:
                continue
                
            if is_safe(board, i, col, n):
                board[i][col] = 1
                
                if solve_util(board, col + 1):
                    return True
                
                board[i][col] = 0
        
        return False
    
    if solve_util(board, 0) == False:
        print("Solution does not exist")
        return None
    
    # Convert solution to required format
    queens = []
    for i in range(n):
        for j in range(n):
            if board[i][j] == 1:
                queens.append(f"{i} {j}")
    
    print(", ".join(queens))

solve_n_queens()