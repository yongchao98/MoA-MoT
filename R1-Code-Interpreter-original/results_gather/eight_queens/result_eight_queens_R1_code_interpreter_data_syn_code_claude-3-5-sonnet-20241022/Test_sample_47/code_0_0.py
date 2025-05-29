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

def solve_n_queens(board, col, n):
    if col >= n:
        return True

    # Skip if there's already a queen in this column
    queen_in_col = False
    for i in range(n):
        if board[i][col] == 1:
            queen_in_col = True
            break
    if queen_in_col:
        return solve_n_queens(board, col + 1, n)

    for i in range(n):
        # Skip if cell is blocked (marked with 2)
        if board[i][col] == 2:
            continue
        
        # Skip if cell already has a queen
        if board[i][col] == 1:
            continue
            
        if is_safe(board, i, col, n):
            board[i][col] = 1
            
            if solve_n_queens(board, col + 1, n):
                return True
                
            board[i][col] = 0
    
    return False

# Initialize the board
n = 8
board = [
    [0,0,0,0,0,0,0,2],  # X marked as 2
    [0,0,0,1,0,0,0,0],
    [0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0],
    [2,0,0,2,0,0,0,0],
    [0,0,2,0,0,0,0,0],
    [0,0,0,0,0,0,0,0]
]

# Solve the puzzle
if solve_n_queens(board, 0, n):
    # Print positions of queens
    positions = []
    for i in range(n):
        for j in range(n):
            if board[i][j] == 1:
                positions.append(f"{i} {j}")
    print(", ".join(positions))
else:
    print("No solution exists")