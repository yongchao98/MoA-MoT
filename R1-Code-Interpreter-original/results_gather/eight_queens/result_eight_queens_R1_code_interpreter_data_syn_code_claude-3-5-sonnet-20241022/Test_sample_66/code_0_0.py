def is_safe(board, row, col, n=8):
    # Check if position is forbidden (X)
    if row == 3 and col == 6:
        return False
    
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
    if col >= n:
        return True
    
    for row in range(n):
        if is_safe(board, row, col):
            board[row][col] = 1
            
            if solve_n_queens(board, col + 1):
                return True
            
            board[row][col] = 0
    
    return False

def get_queen_positions(board, n=8):
    positions = []
    for i in range(n):
        for j in range(n):
            if board[i][j] == 1:
                positions.append(f"{i} {j}")
    return positions

# Initialize board
board = [[0 for x in range(8)] for y in range(8)]

# Solve the puzzle
if solve_n_queens(board, 0):
    positions = get_queen_positions(board)
    print("<<<" + ", ".join(positions) + ">>>")
else:
    print("No solution exists")