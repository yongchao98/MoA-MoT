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
    for i, j in zip(range(row, -1, -1), range(col, -1, -1)):
        if board[i][j] == 1:
            return False
    
    # Check upper right diagonal
    for i, j in zip(range(row, -1, -1), range(col, n)):
        if board[i][j] == 1:
            return False
    
    # Check lower left diagonal
    for i, j in zip(range(row, n), range(col, -1, -1)):
        if board[i][j] == 1:
            return False
    
    # Check lower right diagonal
    for i, j in zip(range(row, n), range(col, n)):
        if board[i][j] == 1:
            return False
    
    return True

def solve_n_queens(n=8):
    # Initialize board with given queen at (2,0)
    board = [[0 for x in range(n)] for y in range(n)]
    board[2][0] = 1  # Given queen
    
    def solve_util(board, col, queens_placed):
        if queens_placed == n:
            return True
            
        # If we're at column 0 and already processed it (due to given queen), skip
        if col == 0:
            return solve_util(board, col + 1, queens_placed)
            
        for i in range(n):
            if is_safe(board, i, col, n):
                board[i][col] = 1
                if solve_util(board, col + 1, queens_placed + 1):
                    return True
                board[i][col] = 0
        
        return False
    
    # Start solving from column 1 since column 0 has the given queen
    if solve_util(board, 1, 1):  # Start with 1 queen already placed
        # Get positions of all queens
        positions = []
        for i in range(n):
            for j in range(n):
                if board[i][j] == 1:
                    positions.append(f"{i} {j}")
        print(", ".join(positions))
    else:
        print("No solution exists")

solve_n_queens()