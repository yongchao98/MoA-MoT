def is_safe(board, row, col, n):
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

def solve_n_queens(board, n):
    positions = []
    
    def solve(row):
        if row >= n:
            return True
        
        # Skip if row already has a queen
        if 1 in board[row]:
            positions.append((row, board[row].index(1)))
            return solve(row + 1)
        
        for col in range(n):
            # Skip if position is forbidden
            if board[row][col] == 'X':
                continue
                
            if is_safe(board, row, col, n):
                board[row][col] = 1
                positions.append((row, col))
                
                if solve(row + 1):
                    return True
                    
                board[row][col] = 0
                positions.pop()
        
        return False
    
    # Initialize board with given state
    board = [
        [0,1,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,'X',0],
        [0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0],
        [0,0,'X','X',0,0,0,0],
        [0,0,0,0,0,0,0,0]
    ]
    
    if solve(0):
        # Sort positions by row for consistent output
        positions.sort(key=lambda x: x[0])
        result = ", ".join(f"{pos[0]} {pos[1]}" for pos in positions)
        print(f"<<{result}>>")
    else:
        print("No solution exists")

solve_n_queens([], 8)