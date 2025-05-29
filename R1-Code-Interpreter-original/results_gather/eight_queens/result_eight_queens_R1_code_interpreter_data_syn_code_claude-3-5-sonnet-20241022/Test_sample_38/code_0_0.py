def is_safe(board, row, col, queens):
    # Check if a queen can be placed on board[row][col]
    
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
    for i, j in zip(range(row, 8, 1), range(col, 8, 1)):
        if board[i][j] == 1:
            return False
    
    # Check diagonal (top-right to bottom-left)
    for i, j in zip(range(row, -1, -1), range(col, 8, 1)):
        if board[i][j] == 1:
            return False
    for i, j in zip(range(row, 8, 1), range(col, -1, -1)):
        if board[i][j] == 1:
            return False
    
    return True

def solve_queens(board, queens_left, queens_positions):
    if queens_left == 0:
        return True, queens_positions
    
    for row in range(8):
        for col in range(8):
            # Skip if cell is not empty or is marked with X
            if board[row][col] != 0 or (row == 2 and col == 0):
                continue
            
            if is_safe(board, row, col, queens_positions):
                board[row][col] = 1
                queens_positions.append((row, col))
                
                success, positions = solve_queens(board, queens_left - 1, queens_positions)
                if success:
                    return True, positions
                
                board[row][col] = 0
                queens_positions.pop()
    
    return False, queens_positions

# Initialize the board
board = [
    [0,0,0,0,0,0,0,0],
    [1,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0],  # X at (2,0)
    [0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0],
    [0,1,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0]
]

# Initialize with existing queens
queens_positions = [(1,0), (5,1)]
success, final_positions = solve_queens(board, 6, queens_positions)

# Format the output
if success:
    result = ", ".join(f"{row} {col}" for row, col in sorted(final_positions))
    print(f"<<{result}>>")
else:
    print("No solution exists")