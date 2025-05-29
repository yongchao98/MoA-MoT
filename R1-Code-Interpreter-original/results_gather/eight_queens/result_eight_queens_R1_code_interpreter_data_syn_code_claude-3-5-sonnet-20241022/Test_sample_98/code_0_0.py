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
    # Initialize board
    board = [[0 for x in range(n)] for y in range(n)]
    
    # Mark restricted positions
    restricted = [(3,2), (7,1)]
    for row, col in restricted:
        board[row][col] = 'X'
    
    def solve_util(board, col):
        if col >= n:
            return True
        
        for i in range(n):
            if board[i][col] != 'X' and is_safe(board, i, col, n):
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
    for i in range(n):
        for j in range(n):
            if board[i][j] == 1:
                positions.append(f"{i} {j}")
    
    return ", ".join(positions)

result = solve_n_queens()
print(result)