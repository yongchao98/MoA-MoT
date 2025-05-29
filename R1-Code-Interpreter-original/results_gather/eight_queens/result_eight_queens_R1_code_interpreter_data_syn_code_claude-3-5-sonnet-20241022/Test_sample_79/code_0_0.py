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

def solve_n_queens(board, col, n, forbidden):
    # Base case: If all queens are placed, return True
    if col >= n:
        return True
    
    # Skip column 7 as it already has a queen
    if col == 7:
        return solve_n_queens(board, col + 1, n, forbidden)
    
    # Try placing queen in each row of the current column
    for row in range(n):
        # Check if current position is forbidden
        if (row, col) in forbidden:
            continue
            
        # Check if queen can be placed
        if is_safe(board, row, col, n):
            board[row][col] = 1
            
            # Recur to place rest of the queens
            if solve_n_queens(board, col + 1, n, forbidden):
                return True
            
            # If placing queen doesn't lead to a solution, backtrack
            board[row][col] = 0
    
    return False

# Initialize the board
n = 8
board = [[0 for x in range(n)] for y in range(n)]
# Place the existing queen
board[0][7] = 1
# Mark forbidden position
forbidden = {(4,1)}

# Solve the puzzle
if solve_n_queens(board, 0, n, forbidden):
    # Collect queen positions
    positions = []
    for i in range(n):
        for j in range(n):
            if board[i][j] == 1:
                positions.append(f"{i} {j}")
    print(", ".join(positions))
else:
    print("No solution exists")