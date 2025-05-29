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
    # Initialize board with constraints
    board = [[0 for x in range(n)] for y in range(n)]
    
    # Mark forbidden positions (X)
    forbidden = [(1,3), (1,5), (3,6), (5,6)]
    for row, col in forbidden:
        board[row][col] = 2  # 2 represents 'X'
    
    def solve_util(board, col):
        if col >= n:
            return True
        
        for i in range(n):
            if board[i][col] != 2 and is_safe(board, i, col, n):  # Check if not forbidden and safe
                board[i][col] = 1
                
                if solve_util(board, col + 1):
                    return True
                
                board[i][col] = 0
        
        return False
    
    if solve_util(board, 0) == False:
        print("Solution does not exist")
        return
    
    # Collect queen positions
    positions = []
    for i in range(n):
        for j in range(n):
            if board[i][j] == 1:
                positions.append(f"{i} {j}")
    
    print("<<<" + ", ".join(positions) + ">>>")

solve_n_queens()