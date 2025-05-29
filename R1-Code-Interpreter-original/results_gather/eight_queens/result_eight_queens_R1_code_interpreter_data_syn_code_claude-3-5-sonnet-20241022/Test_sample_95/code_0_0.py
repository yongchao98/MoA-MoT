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

def solve_n_queens(board, col, n, blocked):
    if col >= n:
        return True
    
    # Try placing queen in each row of current column
    for row in range(n):
        # Check if position is blocked
        if (row, col) in blocked:
            continue
            
        if is_safe(board, row, col, n):
            board[row][col] = 1
            
            if solve_n_queens(board, col + 1, n, blocked):
                return True
                
            board[row][col] = 0
    
    return False

def get_queen_positions(board, n):
    positions = []
    for i in range(n):
        for j in range(n):
            if board[i][j] == 1:
                positions.append(f"{i} {j}")
    return positions

# Initialize board and blocked positions
n = 8
board = [[0 for x in range(n)] for y in range(n)]
blocked = {(0,1), (6,6), (7,6)}  # X positions

# Solve the puzzle
if solve_n_queens(board, 0, n, blocked):
    positions = get_queen_positions(board, n)
    print(", ".join(positions))
else:
    print("No solution exists")