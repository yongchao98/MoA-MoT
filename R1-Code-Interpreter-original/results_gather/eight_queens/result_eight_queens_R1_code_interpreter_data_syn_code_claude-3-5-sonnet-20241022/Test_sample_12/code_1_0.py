def is_safe(board, row, col, n):
    # Check row on both sides
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
            # Check if position is on any diagonal of current position
            if (i + j == row + col) or (i - j == row - col):
                if board[i][j] == 1:
                    return False
    
    return True

def solve_n_queens(n=8):
    # Initialize board with pre-placed queen and X
    board = [[0 for x in range(n)] for y in range(n)]
    board[2][7] = 1  # Pre-placed queen
    board[2][5] = -1  # X position
    
    def solve_util(board, queens_placed):
        if queens_placed == n:
            return True
            
        # If we're placing the first queen, skip the pre-placed queen's count
        if queens_placed == 0:
            queens_placed = 1
            
        for col in range(n):
            for row in range(n):
                # Skip if position is X or the pre-placed queen
                if (row == 2 and col == 5) or (row == 2 and col == 7):
                    continue
                    
                if board[row][col] == 0 and is_safe(board, row, col, n):
                    board[row][col] = 1
                    
                    if solve_util(board, queens_placed + 1):
                        return True
                    
                    board[row][col] = 0
        
        return False
    
    if solve_util(board, 0) == False:
        print("Solution does not exist")
        return None
    
    # Collect queen positions
    queens = []
    for i in range(n):
        for j in range(n):
            if board[i][j] == 1:
                queens.append(f"{i} {j}")
    
    print("<<<" + ", ".join(queens) + ">>>")

solve_n_queens()