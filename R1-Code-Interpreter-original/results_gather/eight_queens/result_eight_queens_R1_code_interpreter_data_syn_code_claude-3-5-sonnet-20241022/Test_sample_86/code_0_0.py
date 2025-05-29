def is_safe(board, row, col, n):
    # Check row on left side
    for j in range(col):
        if board[row][j] == 1:
            return False
    
    # Check upper diagonal on left side
    for i, j in zip(range(row-1, -1, -1), range(col-1, -1, -1)):
        if board[i][j] == 1:
            return False
    
    # Check lower diagonal on left side
    for i, j in zip(range(row+1, n), range(col-1, -1, -1)):
        if board[i][j] == 1:
            return False
    
    return True

def solve_queens(board, col, n, blocked_pos):
    # Base case: If all queens are placed, return True
    if col >= n:
        return True
    
    # Consider this column and try placing queens in all rows one by one
    for row in range(n):
        # Check if we can place a queen here
        if (row, col) != blocked_pos and board[row][col] != 1 and is_safe(board, row, col, n):
            # Place the queen
            board[row][col] = 1
            
            # Recur to place rest of the queens
            if solve_queens(board, col + 1, n, blocked_pos):
                return True
            
            # If placing queen doesn't lead to a solution, backtrack
            board[row][col] = 0
    
    return False

def find_queens_positions():
    n = 8
    # Initialize board with given state
    board = [
        [0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,1,0],
        [0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0],
        [0,0,0,-1,0,0,0,0],  # -1 represents blocked position
        [0,0,0,0,0,0,0,0]
    ]
    
    blocked_pos = (6, 3)
    
    # Start from column 0
    if solve_queens(board, 0, n, blocked_pos):
        # Collect positions of queens
        positions = []
        for i in range(n):
            for j in range(n):
                if board[i][j] == 1:
                    positions.append(f"{i} {j}")
        print("<<<" + ", ".join(positions) + ">>>")
    else:
        print("No solution exists")

find_queens_positions()