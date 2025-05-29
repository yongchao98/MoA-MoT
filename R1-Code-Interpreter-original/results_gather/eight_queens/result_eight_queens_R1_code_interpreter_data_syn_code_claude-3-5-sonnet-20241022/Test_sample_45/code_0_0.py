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
    # Initialize board with existing queen and X positions
    board = [[0 for x in range(n)] for y in range(n)]
    
    # Mark X positions
    x_positions = [(1,6), (3,5), (6,4), (6,6)]
    for x, y in x_positions:
        board[x][y] = 2  # 2 represents 'X'
    
    # Place existing queen
    board[3][7] = 1
    
    def solve_util(board, col):
        if col >= n:
            return True
        
        # Skip column 7 as it already has a queen
        if col == 7:
            return solve_util(board, col + 1)
        
        for i in range(n):
            if board[i][col] != 2 and is_safe(board, i, col, n):  # Check if not 'X' and safe
                board[i][col] = 1
                
                if solve_util(board, col + 1):
                    return True
                
                board[i][col] = 0
        
        return False
    
    if solve_util(board, 0):
        # Get queen positions
        positions = []
        for i in range(n):
            for j in range(n):
                if board[i][j] == 1:
                    positions.append(f"{i} {j}")
        print("<<<" + ", ".join(positions) + ">>>")
    else:
        print("No solution exists")

solve_n_queens()