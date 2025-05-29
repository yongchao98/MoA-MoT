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

def solve_queens(board, queens_placed, queens):
    if queens_placed == 8:
        return True
    
    # Try placing queen in each row
    for row in range(8):
        # Try placing queen in each column
        for col in range(8):
            # Skip if queen is already placed here
            if board[row][col] == 1:
                continue
            
            # Check if we can place a queen here
            if is_safe(board, row, col, queens):
                # Place the queen
                board[row][col] = 1
                queens.append((row, col))
                
                # Recursively place remaining queens
                if solve_queens(board, queens_placed + 1, queens):
                    return True
                
                # If placing queen here didn't work, remove it
                board[row][col] = 0
                queens.pop()
    
    return False

# Initialize board with the given queen
board = [[0 for _ in range(8)] for _ in range(8)]
board[7][7] = 1  # Pre-placed queen
queens = [(7, 7)]  # Keep track of queen positions

# Solve the puzzle
if solve_queens(board, 1, queens):
    # Sort queens by row for consistent output
    queens.sort()
    # Format the output
    result = ", ".join(f"{row} {col}" for row, col in queens)
    print(f"<<{result}>>")
else:
    print("No solution exists")