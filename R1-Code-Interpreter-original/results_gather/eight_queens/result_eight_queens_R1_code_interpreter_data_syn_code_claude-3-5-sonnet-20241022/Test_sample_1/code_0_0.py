def is_safe(board, row, col, n):
    # Check row
    for j in range(n):
        if board[row][j] == 1:
            return False
    
    # Check column
    for i in range(n):
        if board[i][col] == 1:
            return False
    
    # Check upper left diagonal
    for i, j in zip(range(row-1, -1, -1), range(col-1, -1, -1)):
        if board[i][j] == 1:
            return False
    
    # Check upper right diagonal
    for i, j in zip(range(row-1, -1, -1), range(col+1, n)):
        if board[i][j] == 1:
            return False
    
    # Check lower left diagonal
    for i, j in zip(range(row+1, n), range(col-1, -1, -1)):
        if board[i][j] == 1:
            return False
    
    # Check lower right diagonal
    for i, j in zip(range(row+1, n), range(col+1, n)):
        if board[i][j] == 1:
            return False
    
    return True

def solve_n_queens(n=8):
    # Initialize board with existing queen
    board = [[0 for x in range(n)] for y in range(n)]
    board[2][4] = 1  # Given queen position
    queens = [(2,4)]  # Store queen positions
    
    def solve_util(board, queens_placed):
        if queens_placed == n:
            return True
        
        for row in range(n):
            for col in range(n):
                if board[row][col] == 0 and is_safe(board, row, col, n):
                    board[row][col] = 1
                    queens.append((row,col))
                    if solve_util(board, queens_placed + 1):
                        return True
                    board[row][col] = 0
                    queens.pop()
        return False
    
    if solve_util(board, 1):  # Start with 1 as we already have one queen
        # Sort queens by row for consistent output
        queens.sort()
        result = ", ".join(f"{row} {col}" for row, col in queens)
        print(f"<<{result}>>")
    else:
        print("No solution exists")

solve_n_queens()