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

def solve_queens(board, queens_placed=0, queens_positions=None):
    if queens_positions is None:
        queens_positions = []
        # Add initial queens
        for i in range(8):
            for j in range(8):
                if board[i][j] == 1:
                    queens_positions.append((i, j))
                    queens_placed += 1
    
    if queens_placed == 8:
        # Sort positions by row for consistent output
        queens_positions.sort()
        result = ", ".join(f"{row} {col}" for row, col in queens_positions)
        print(f"<<<{result}>>>")
        return True
    
    for row in range(8):
        for col in range(8):
            if board[row][col] == 0 and is_safe(board, row, col, queens_positions):
                board[row][col] = 1
                queens_positions.append((row, col))
                
                if solve_queens(board, queens_placed + 1, queens_positions):
                    return True
                
                board[row][col] = 0
                queens_positions.pop()
    
    return False

# Initial board
board = [
    [0,0,0,0,0,1,0,0],
    [0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0],
    [0,0,0,1,0,0,0,0],
    [0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0]
]

solve_queens(board)