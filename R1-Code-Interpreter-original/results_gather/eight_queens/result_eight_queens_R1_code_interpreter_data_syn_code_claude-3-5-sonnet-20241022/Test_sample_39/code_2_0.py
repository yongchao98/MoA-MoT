def is_safe(board, row, col, n):
    # Check row
    for j in range(n):
        if board[row][j] == 1:
            return False
    
    # Check column
    for i in range(n):
        if board[i][col] == 1:
            return False
    
    # Check diagonals
    for i in range(n):
        for j in range(n):
            if board[i][j] == 1:
                if abs(row - i) == abs(col - j):
                    return False
    
    return True

def solve_n_queens(n=8):
    # Initialize board with given constraints
    board = [
        [0,0,0,0,0,1,'X',0],
        [0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0],
        ['X',0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0],
        ['X',0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0],
        [0,0,0,'X',0,0,0,0]
    ]
    
    def solve_util(board, col):
        # If all queens are placed, return True
        if col >= n:
            return True
            
        # Skip column 5 as queen is already placed there
        if col == 5:
            return solve_util(board, col + 1)
            
        for i in range(n):
            # Skip if position has 'X'
            if board[i][col] == 'X':
                continue
                
            if is_safe(board, i, col, n):
                # Place queen
                board[i][col] = 1
                
                # Recur to place rest of the queens
                if solve_util(board, col + 1):
                    return True
                    
                # If placing queen doesn't lead to a solution, remove queen
                board[i][col] = 0
        
        return False
    
    if solve_util(board, 0):
        # Collect queen positions
        positions = []
        for i in range(n):
            for j in range(n):
                if board[i][j] == 1:
                    positions.append(f"{i} {j}")
        print(",".join(sorted(positions)))
    else:
        print("No solution exists")

solve_n_queens()