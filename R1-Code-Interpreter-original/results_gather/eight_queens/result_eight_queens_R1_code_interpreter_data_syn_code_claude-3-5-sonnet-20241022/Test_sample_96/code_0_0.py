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
    # Initialize board with forbidden positions marked as 2
    board = [[0 for x in range(n)] for y in range(n)]
    # Mark forbidden positions (X)
    forbidden = [(2,6), (3,1), (3,6), (6,1)]
    for x, y in forbidden:
        board[x][y] = 2
    
    def solve_util(board, col):
        if col >= n:
            return True
        
        for i in range(n):
            if board[i][col] != 2 and is_safe(board, i, col, n):
                board[i][col] = 1
                
                if solve_util(board, col + 1):
                    return True
                
                board[i][col] = 0
        
        return False
    
    if solve_util(board, 0) == False:
        print("Solution does not exist")
        return None
    
    # Get queen positions
    positions = []
    for col in range(n):
        for row in range(n):
            if board[row][col] == 1:
                positions.append(f"{row} {col}")
    
    return positions

# Solve and print result
result = solve_n_queens()
if result:
    print(", ".join(result))