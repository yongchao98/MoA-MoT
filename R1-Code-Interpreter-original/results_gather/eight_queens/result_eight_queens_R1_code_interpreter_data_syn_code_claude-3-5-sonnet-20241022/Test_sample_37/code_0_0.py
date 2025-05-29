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
    
    # Check if position is marked with X (2)
    if board[row][col] == 2:
        return False
        
    return True

def solve_n_queens(board, col, n):
    if col >= n:
        return True

    # Skip column if it already has a queen
    for i in range(n):
        if board[i][col] == 1:
            return solve_n_queens(board, col + 1, n)

    for i in range(n):
        if is_safe(board, i, col, n):
            # Place queen
            board[i][col] = 1
            
            if solve_n_queens(board, col + 1, n):
                return True
                
            # If placing queen doesn't lead to a solution, remove it
            board[i][col] = 0
    
    return False

# Initialize board with given configuration
board = [
    [0,0,0,0,0,0,0,0],
    [0,2,0,0,0,0,0,0],
    [0,0,0,0,2,0,0,0],
    [0,0,0,0,0,0,1,0],
    [0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0],
    [0,0,2,0,0,0,0,0]
]

# Solve the puzzle
if solve_n_queens(board, 0, 8):
    # Get queen positions
    positions = []
    for i in range(8):
        for j in range(8):
            if board[i][j] == 1:
                positions.append(f"{i} {j}")
    print(", ".join(positions))
else:
    print("No solution exists")