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
    
    # Consider this column and try placing this queen in all rows one by one
    for i in range(n):
        # Check if position is forbidden
        if (i, col) in forbidden:
            continue
            
        # Check if queen can be placed on board[i][col]
        if is_safe(board, i, col, n):
            # Place this queen in board[i][col]
            board[i][col] = 1
            
            # Recur to place rest of the queens
            if solve_n_queens(board, col + 1, n, forbidden):
                return True
            
            # If placing queen in board[i][col] doesn't lead to a solution,
            # then remove queen from board[i][col]
            board[i][col] = 0
    
    # If queen can't be placed in any row in this column col, return False
    return False

def get_queen_positions(board, n):
    positions = []
    for i in range(n):
        for j in range(n):
            if board[i][j] == 1:
                positions.append(f"{i} {j}")
    return positions

# Initialize the board
n = 8
board = [
    [0,0,0,0,0,0,0,0],
    ['X',0,0,0,0,0,'X',0],
    [0,0,0,0,'X',0,0,0],
    [0,0,0,0,0,0,0,0],
    [0,0,0,'X',0,0,0,0],
    [0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0]
]

# Convert board to numeric (0s and 1s)
numeric_board = [[0 if cell == 0 else 0 for cell in row] for row in board]

# Create set of forbidden positions
forbidden = set()
for i in range(n):
    for j in range(n):
        if board[i][j] == 'X':
            forbidden.add((i, j))

# Solve the puzzle
if solve_n_queens(numeric_board, 0, n, forbidden):
    positions = get_queen_positions(numeric_board, n)
    print("<<<" + ", ".join(positions) + ">>>")
else:
    print("No solution exists")