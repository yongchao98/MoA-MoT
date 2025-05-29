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
    
    # Check if position is blocked
    if board[row][col] == 'X':
        return False
        
    return True

def solve_n_queens(board, col, n):
    # Base case: If all queens are placed, return True
    if col >= n:
        return True
    
    # Skip column 3 as it already has a queen in row 0
    if col == 3:
        return solve_n_queens(board, col + 1, n)
    
    # Consider this column and try placing this queen in all rows one by one
    for i in range(n):
        if is_safe(board, i, col, n):
            # Place this queen in board[i][col]
            board[i][col] = 1
            
            # Recur to place rest of the queens
            if solve_n_queens(board, col + 1, n):
                return True
            
            # If placing queen in board[i][col] doesn't lead to a solution,
            # then remove queen from board[i][col]
            board[i][col] = 0
    
    # If queen can't be placed in any row in this column col, return False
    return False

# Initialize the board
n = 8
board = [[0 for x in range(n)] for y in range(n)]

# Place the pre-existing queen
board[0][3] = 1

# Place the blocked position
board[1][5] = 'X'

# Solve the puzzle
if solve_n_queens(board, 0, n):
    # Collect queen positions
    positions = []
    for i in range(n):
        for j in range(n):
            if board[i][j] == 1:
                positions.append(f"{i} {j}")
    print(", ".join(positions))
else:
    print("No solution exists")