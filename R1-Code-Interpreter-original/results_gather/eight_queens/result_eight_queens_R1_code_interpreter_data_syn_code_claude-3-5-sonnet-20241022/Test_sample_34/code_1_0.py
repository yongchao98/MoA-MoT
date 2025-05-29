def is_safe(board, row, col, n):
    # Check if position has 'X'
    if board[row][col] == 'X':
        return False
        
    # Check row
    for j in range(n):
        if board[row][j] == 1:
            return False
    
    # Check column
    for i in range(n):
        if board[i][col] == 1:
            return False
    
    # Check diagonal (top-left to bottom-right)
    for i, j in zip(range(row, -1, -1), range(col, -1, -1)):
        if board[i][j] == 1:
            return False
    for i, j in zip(range(row, n, 1), range(col, n, 1)):
        if board[i][j] == 1:
            return False
    
    # Check diagonal (top-right to bottom-left)
    for i, j in zip(range(row, -1, -1), range(col, n, 1)):
        if board[i][j] == 1:
            return False
    for i, j in zip(range(row, n, 1), range(col, -1, -1)):
        if board[i][j] == 1:
            return False
    
    return True

def solve_n_queens(n=8):
    # Initialize board with given configuration
    board = [
        [0,0,0,0,1,0,0,0],
        [0,1,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0],
        [0,0,'X',0,0,0,0,0],
        [0,0,'X',0,0,0,0,0],
        [0,0,0,0,0,0,0,0]
    ]
    
    def count_queens():
        return sum(row.count(1) for row in board)
    
    def solve_util(row):
        if count_queens() == 8:
            return True
            
        if row >= n:
            return False
            
        for col in range(n):
            if board[row][col] != 1 and is_safe(board, row, col, n):
                board[row][col] = 1
                if solve_util(row + 1):
                    return True
                board[row][col] = 0
        
        return solve_util(row + 1) if any(board[row][j] == 1 for j in range(n)) else False
    
    if solve_util(0):
        # Get positions of queens
        positions = []
        for i in range(n):
            for j in range(n):
                if board[i][j] == 1:
                    positions.append(f"{i} {j}")
        print(", ".join(positions))
    else:
        print("No solution exists")

solve_n_queens()