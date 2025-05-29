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

def solve_queens(board, queens, positions):
    if queens == 8:
        return True, positions
    
    for row in range(8):
        for col in range(8):
            # Skip if cell is forbidden (X) or already has a queen
            if board[row][col] == 'X' or board[row][col] == 1:
                continue
                
            if is_safe(board, row, col, queens):
                board[row][col] = 1
                positions.append(f"{row} {col}")
                success, result = solve_queens(board, queens + 1, positions)
                if success:
                    return True, result
                board[row][col] = 0
                positions.pop()
    
    return False, positions

# Initialize the board
board = [
    [0,0,0,0,0,0,0,0],
    ['X',0,'X',0,0,0,'X',0],
    [0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0],
    [0,0,'X',0,0,0,0,0]
]

success, positions = solve_queens(board, 0, [])
if success:
    print("<<<" + ", ".join(positions) + ">>>")
else:
    print("No solution exists")