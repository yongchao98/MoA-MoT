def is_safe(board, row, col, n=8):
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

def solve_n_queens(board, col, n=8):
    # Base case: If all queens are placed, return True
    if col >= n:
        return True
    
    # Check if current column already has a queen
    queen_in_col = False
    for i in range(n):
        if board[i][col] == 1:
            queen_in_col = True
            # Move to next column
            if solve_n_queens(board, col + 1):
                return True
            break
    
    if queen_in_col:
        return False
    
    # Try placing queen in each row of current column
    for i in range(n):
        # Check if position is marked with X
        if board[i][col] == 'X':
            continue
            
        # Check if queen can be placed
        if is_safe(board, i, col):
            # Place the queen
            board[i][col] = 1
            
            # Recur to place rest of the queens
            if solve_n_queens(board, col + 1):
                return True
            
            # If placing queen doesn't lead to solution, backtrack
            board[i][col] = 0
    
    return False

# Initialize the board with given constraints
board = [
    ['X', 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 1, 0],
    [0, 0, 'X', 0, 0, 0, 0, 0],
    [0, 'X', 1, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0]
]

# Solve the puzzle
if solve_n_queens(board, 0):
    # Get queen positions
    positions = []
    for i in range(8):
        for j in range(8):
            if board[i][j] == 1:
                positions.append(f"{i} {j}")
    print(", ".join(positions))
else:
    print("No solution exists")