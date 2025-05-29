def is_safe(board, row, col):
    # Check row on left side
    for j in range(col):
        if board[row][j] == 1:
            return False
    
    # Check upper diagonal on left side
    i, j = row - 1, col - 1
    while i >= 0 and j >= 0:
        if board[i][j] == 1:
            return False
        i -= 1
        j -= 1
    
    # Check lower diagonal on left side
    i, j = row + 1, col - 1
    while i < 8 and j >= 0:
        if board[i][j] == 1:
            return False
        i += 1
        j -= 1
        
    # Check if this position conflicts with the fixed queen at (7,5)
    if col != 5:  # Skip checking when we're placing in column 5
        # Check same row as fixed queen
        if row == 7:
            return False
        # Check diagonal with fixed queen
        if abs(row - 7) == abs(col - 5):
            return False
    
    return True

def solve_n_queens():
    board = [[0 for _ in range(8)] for _ in range(8)]
    # Place the fixed queen
    board[7][5] = 1
    queens = []
    
    def solve_util(col):
        # If we've placed queens in all columns except 5
        if len(queens) == 7:
            queens.append((7, 5))  # Add the fixed queen
            return True
            
        # Skip column 5 as it already has a queen
        if col == 5:
            return solve_util(col + 1)
            
        # If we've placed all queens
        if col >= 8:
            return True
            
        # Try placing queen in each row of current column
        for row in range(8):
            if is_safe(board, row, col):
                board[row][col] = 1
                queens.append((row, col))
                
                if solve_util(col + 1):
                    return True
                    
                board[row][col] = 0
                queens.pop()
                
        return False
    
    if solve_util(0):
        # Sort queens by row for consistent output
        queens.sort()
        result = ", ".join(f"{row} {col}" for row, col in queens)
        print(f"<<{result}>>")
    else:
        print("No solution exists")

solve_n_queens()