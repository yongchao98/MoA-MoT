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
    # Initialize board with pre-placed queen and X
    board = [[0 for x in range(n)] for y in range(n)]
    board[2][7] = 1  # Pre-placed queen
    board[2][5] = -1  # X position
    
    def solve_util(board, col):
        if col >= n:
            return True
        
        # Skip column 7 as it already has a queen
        if col == 7:
            return solve_util(board, col + 1)
        
        for i in range(n):
            # Skip row 2, col 5 as it's marked with X
            if col == 5 and i == 2:
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
    
    # Collect queen positions
    queens = []
    for i in range(n):
        for j in range(n):
            if board[i][j] == 1:
                queens.append(f"{i} {j}")
    
    print("<<<" + ", ".join(queens) + ">>>")

solve_n_queens()