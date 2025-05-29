def is_safe(board, row, col, forbidden):
    # Check if position is forbidden
    if (row, col) in forbidden:
        return False
    
    # Check row
    for j in range(8):
        if board[row][j] == 1:
            return False
    
    # Check column
    for i in range(8):
        if board[i][col] == 1:
            return False
    
    # Check diagonal (top-left to bottom-right)
    for i, j in zip(range(row-1, -1, -1), range(col-1, -1, -1)):
        if board[i][j] == 1:
            return False
    for i, j in zip(range(row+1, 8), range(col+1, 8)):
        if board[i][j] == 1:
            return False
    
    # Check diagonal (top-right to bottom-left)
    for i, j in zip(range(row-1, -1, -1), range(col+1, 8)):
        if board[i][j] == 1:
            return False
    for i, j in zip(range(row+1, 8), range(col-1, -1, -1)):
        if board[i][j] == 1:
            return False
    
    return True

def solve_queens(board, col, forbidden):
    # Base case: all queens are placed
    if col >= 8:
        return True
    
    # Try placing queen in each row of current column
    for row in range(8):
        if is_safe(board, row, col, forbidden):
            board[row][col] = 1
            
            if solve_queens(board, col + 1, forbidden):
                return True
            
            board[row][col] = 0
    
    return False

# Initialize board and forbidden positions
board = [[0 for _ in range(8)] for _ in range(8)]
forbidden = {(5,1), (5,5), (7,4)}

# Solve the puzzle
if solve_queens(board, 0, forbidden):
    # Get queen positions
    positions = []
    for col in range(8):
        for row in range(8):
            if board[row][col] == 1:
                positions.append(f"{row} {col}")
    print(", ".join(positions))
else:
    print("No solution exists")