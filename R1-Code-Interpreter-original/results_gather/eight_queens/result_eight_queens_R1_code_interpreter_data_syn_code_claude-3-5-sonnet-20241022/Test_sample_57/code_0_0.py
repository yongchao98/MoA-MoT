def is_safe(board, row, col, n):
    # Check if position (5,3) is being considered
    if row == 5 and col == 3:
        return False
        
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

def solve_n_queens(n):
    board = [[0 for x in range(n)] for y in range(n)]
    
    def solve_util(board, col):
        if col >= n:
            return True
        
        for i in range(n):
            if is_safe(board, i, col, n):
                board[i][col] = 1
                
                if solve_util(board, col + 1):
                    return True
                
                board[i][col] = 0
        
        return False
    
    if solve_util(board, 0) == False:
        print("Solution does not exist")
        return False
    
    # Collect queen positions
    positions = []
    for col in range(n):
        for row in range(n):
            if board[row][col] == 1:
                positions.append(f"{row} {col}")
    
    print(", ".join(positions))
    return True

solve_n_queens(8)