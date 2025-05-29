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
    if col >= n:
        return True
    
    for row in range(n):
        if (row, col) in forbidden:
            continue
            
        if is_safe(board, row, col, n):
            board[row][col] = 1
            
            if solve_n_queens(board, col + 1, n, forbidden):
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

# Initialize the board
n = 8
board = [[0 for x in range(n)] for y in range(n)]

# Define forbidden positions (X)
forbidden = {(5,0), (5,1), (5,2)}

# Solve the puzzle
if solve_n_queens(board, 0, n, forbidden):
    positions = get_queen_positions(board, n)
    print("<<<" + ", ".join(positions) + ">>>")
else:
    print("No solution exists")