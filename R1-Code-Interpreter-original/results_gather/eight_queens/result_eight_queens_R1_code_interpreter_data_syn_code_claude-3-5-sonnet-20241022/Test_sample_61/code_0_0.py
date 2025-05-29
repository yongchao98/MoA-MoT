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

def solve_n_queens(board, col, n, restricted):
    if col >= n:
        return True
    
    # Try placing queen in each row of current column
    for row in range(n):
        # Check if position is restricted
        if (row, col) in restricted:
            continue
            
        if is_safe(board, row, col, n):
            board[row][col] = 1
            
            if solve_n_queens(board, col + 1, n, restricted):
                return True
                
            board[row][col] = 0
    
    return False

def get_queen_positions(board):
    positions = []
    for i in range(len(board)):
        for j in range(len(board)):
            if board[i][j] == 1:
                positions.append(f"{i} {j}")
    return positions

# Initialize board and restricted positions
n = 8
board = [[0 for x in range(n)] for y in range(n)]
restricted = [(3,4), (4,6), (6,1)]  # X positions

# Solve the puzzle
if solve_n_queens(board, 0, n, restricted):
    positions = get_queen_positions(board)
    print("<<<" + ", ".join(positions) + ">>>")
else:
    print("No solution exists")