def is_safe(board, row, col, n):
    # Check if a queen can be placed on board[row][col]
    
    # Check this column
    for i in range(n):
        if board[i][col] == 1:
            return False
    
    # Check this row
    for j in range(n):
        if board[row][j] == 1:
            return False
    
    # Check diagonals
    for i in range(n):
        for j in range(n):
            if board[i][j] == 1:
                if abs(row - i) == abs(col - j):
                    return False
    
    return True

def solve_n_queens(n=8):
    # Initialize board with given configuration
    board = [[0 for x in range(n)] for y in range(n)]
    
    # Set initial queen and restricted positions
    board[1][4] = 1  # Initial queen
    restricted = {(4,7), (5,1), (7,0)}  # X positions
    
    queens = [(1,4)]  # Keep track of queen positions
    
    def solve_util(board, queens_placed):
        if queens_placed == 8:
            return True
            
        # Try all rows and columns
        for row in range(n):
            for col in range(n):
                # Skip if position is restricted or already has a queen
                if (row,col) in restricted or board[row][col] == 1:
                    continue
                # Skip if this is the initial queen's position
                if (row == 1 and col == 4):
                    continue
                    
                if is_safe(board, row, col, n):
                    board[row][col] = 1
                    queens.append((row,col))
                    
                    if solve_util(board, queens_placed + 1):
                        return True
                        
                    board[row][col] = 0
                    queens.pop()
        return False
    
    # Start solving with 1 queen already placed
    solve_util(board, 1)
    
    # Format the output
    result = ", ".join(f"{row} {col}" for row, col in sorted(queens))
    print(f"<<<{result}>>>")

solve_n_queens()