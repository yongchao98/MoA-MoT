def is_safe(board, row, col, n):
    # Check row
    for j in range(n):
        if board[row][j] == 1:
            return False
    
    # Check column
    for i in range(n):
        if board[i][col] == 1:
            return False
    
    # Check upper right diagonal
    for i, j in zip(range(row-1, -1, -1), range(col+1, n)):
        if board[i][j] == 1:
            return False
    
    # Check upper left diagonal
    for i, j in zip(range(row-1, -1, -1), range(col-1, -1, -1)):
        if board[i][j] == 1:
            return False
            
    # Check lower right diagonal
    for i, j in zip(range(row+1, n), range(col+1, n)):
        if board[i][j] == 1:
            return False
            
    # Check lower left diagonal
    for i, j in zip(range(row+1, n), range(col-1, -1, -1)):
        if board[i][j] == 1:
            return False
    
    # Check if position is marked with X (2)
    if board[row][col] == 2:
        return False
        
    return True

def solve_n_queens(board, queens_placed, n):
    if queens_placed == n:
        return True

    # Try placing queen in each row and column
    for i in range(n):
        for j in range(n):
            # Skip if there's already a queen or if it's not safe
            if board[i][j] != 1 and is_safe(board, i, j, n):
                # Place queen
                temp = board[i][j]
                board[i][j] = 1
                
                if solve_n_queens(board, queens_placed + 1, n):
                    return True
                    
                # If placing queen doesn't lead to a solution, remove it
                board[i][j] = temp
    
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

# Count existing queens
initial_queens = sum(row.count(1) for row in board)

# Solve the puzzle
if solve_n_queens(board, initial_queens, 8):
    # Get queen positions
    positions = []
    for i in range(8):
        for j in range(8):
            if board[i][j] == 1:
                positions.append(f"{i} {j}")
    print(", ".join(sorted(positions)))
else:
    print("No solution exists")