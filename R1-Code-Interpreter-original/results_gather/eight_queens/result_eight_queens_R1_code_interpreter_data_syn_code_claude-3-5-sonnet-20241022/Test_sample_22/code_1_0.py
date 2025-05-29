def is_safe(board, row, col):
    n = 8
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
            # If position has a queen and is on the same diagonal
            if board[i][j] == 1:
                if abs(row - i) == abs(col - j):
                    return False
    
    return True

def solve_queens():
    n = 8
    # Initialize board with pre-placed queen and X
    board = [[0 for _ in range(n)] for _ in range(n)]
    board[1][4] = 1  # Pre-placed queen
    board[3][0] = 'X'  # Forbidden position
    queens_placed = 1  # Start with 1 as we have one queen pre-placed
    
    def solve_util(board, queens_placed):
        if queens_placed == 8:
            return True
            
        for i in range(n):
            for j in range(n):
                # Skip if cell already has queen or is forbidden
                if board[i][j] != 0:
                    continue
                    
                # Skip the pre-placed queen's position
                if i == 1 and j == 4:
                    continue
                    
                # Try placing a queen
                if is_safe(board, i, j):
                    board[i][j] = 1
                    if solve_util(board, queens_placed + 1):
                        return True
                    board[i][j] = 0
        
        return False
    
    if solve_util(board, queens_placed):
        # Collect and print queen positions
        positions = []
        for i in range(n):
            for j in range(n):
                if board[i][j] == 1:
                    positions.append(f"{i} {j}")
        print(", ".join(positions))
    else:
        print("No solution exists")

solve_queens()