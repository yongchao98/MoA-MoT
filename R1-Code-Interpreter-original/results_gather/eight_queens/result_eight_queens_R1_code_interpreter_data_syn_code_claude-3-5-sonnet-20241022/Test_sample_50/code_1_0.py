def is_safe(board, row, col, n):
    # Check row on left side
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
    # Initialize board with the given queen at (6,2)
    board = [[0 for x in range(n)] for y in range(n)]
    board[6][2] = 1  # Given queen
    queens = [(6,2)]  # Keep track of placed queens
    
    def solve_util(row=0):
        if len(queens) == n:
            return True
        
        # Skip row 6 as it already has a queen
        if row == 6:
            return solve_util(row + 1)
        
        if row >= n:
            return False
        
        # Try placing queen in current row
        for col in range(n):
            if is_safe(board, row, col, n):
                board[row][col] = 1
                queens.append((row, col))
                
                if solve_util(row + 1):
                    return True
                    
                board[row][col] = 0
                queens.pop()
                
        return False
    
    if solve_util():
        # Sort queens by row for consistent output
        queens.sort()
        print(",".join(f"{row} {col}" for row, col in queens))
    else:
        print("No solution exists")

solve_n_queens()