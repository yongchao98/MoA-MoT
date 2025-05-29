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

def solve_n_queens(n=8):
    # Initialize board with given queens and blocked position
    board = [
        [0,1,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0],
        [0,0,1,0,0,0,0,0],
        [0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,-1],  # -1 represents blocked position
        [0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0]
    ]
    
    def solve_util(board, queens_placed):
        if queens_placed == 8:
            return True
            
        for i in range(n):
            for j in range(n):
                if board[i][j] == 0 and is_safe(board, i, j, n):
                    board[i][j] = 1
                    if solve_util(board, queens_placed + 1):
                        return True
                    board[i][j] = 0
        return False
    
    # Count initial queens
    initial_queens = sum(row.count(1) for row in board)
    
    if solve_util(board, initial_queens):
        # Collect queen positions
        positions = []
        for i in range(n):
            for j in range(n):
                if board[i][j] == 1:
                    positions.append(f"{i} {j}")
        print(", ".join(positions))
    else:
        print("No solution exists")

solve_n_queens()