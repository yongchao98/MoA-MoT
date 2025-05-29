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

def solve_queens(board, col, queens, forbidden, result):
    # Base case: If all queens are placed, return True
    if len(queens) == 8:
        result.extend(queens)
        return True
    
    # Try placing queen in each row of the current column
    for row in range(8):
        if (row, col) not in forbidden and is_safe(board, row, col, queens):
            board[row][col] = 1
            queens.append((row, col))
            
            # Recur to place rest of the queens
            if solve_queens(board, col + 1 if col < 7 else 0, queens, forbidden, result):
                return True
            
            # If placing queen at (row, col) didn't lead to a solution,
            # remove queen from board and try next position
            board[row][col] = 0
            queens.pop()
    
    return False

# Initialize the board and forbidden positions
board = [[0 for _ in range(8)] for _ in range(8)]
forbidden = {(2,1), (3,7), (5,7), (6,1)}
queens = []
result = []

# Try placing queens starting from first column
solve_queens(board, 0, queens, forbidden, result)

# Format and print the result
formatted_result = ", ".join(f"{row} {col}" for row, col in result)
print(f"<<<{formatted_result}>>>")