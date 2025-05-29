def is_safe(board, row, col, queens):
    # Check row
    for j in range(8):
        if board[row][j] == 1:
            return False
    
    # Check column
    for i in range(8):
        if board[i][col] == 1:
            return False
    
    # Check diagonal (top-left to bottom-right)
    for i, j in zip(range(row, -1, -1), range(col, -1, -1)):
        if board[i][j] == 1:
            return False
    for i, j in zip(range(row, 8), range(col, 8)):
        if board[i][j] == 1:
            return False
    
    # Check diagonal (top-right to bottom-left)
    for i, j in zip(range(row, -1, -1), range(col, 8)):
        if board[i][j] == 1:
            return False
    for i, j in zip(range(row, 8), range(col, -1, -1)):
        if board[i][j] == 1:
            return False
    
    return True

def solve_queens(board, queens, solutions):
    if queens == []:
        # Initialize empty list of queen positions
        queens = []
    
    if len(queens) == 8:
        # Found a solution
        queens.sort()  # Sort by row number
        result = ", ".join(f"{row} {col}" for row, col in queens)
        print(f"<<<{result}>>>")
        return True
    
    for row in range(8):
        for col in range(8):
            # Skip if position is forbidden (X) or already has a queen
            if board[row][col] == 'X' or board[row][col] == 1:
                continue
            
            # Skip if this position is already in queens list
            if (row, col) in queens:
                continue
            
            if is_safe(board, row, col, queens):
                board[row][col] = 1
                queens.append((row, col))
                
                if solve_queens(board, queens, solutions):
                    return True
                
                board[row][col] = 0
                queens.pop()
    
    return False

# Initialize the board
board = [
    [0,0,0,0,0,0,0,0],
    ['X',0,0,0,'X',0,0,0],
    [0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0],
    [0,0,0,'X',0,0,0,0],
    [0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0],
    [0,0,'X',0,0,0,0,0]
]

solutions = []
solve_queens(board, [], solutions)