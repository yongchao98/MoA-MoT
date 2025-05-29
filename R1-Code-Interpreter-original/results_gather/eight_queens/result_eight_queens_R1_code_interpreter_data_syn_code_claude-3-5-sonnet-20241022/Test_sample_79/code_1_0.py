def is_safe(board, row, col, n):
    # Check row on both sides
    for j in range(n):
        if board[row][j] == 1:
            return False
    
    # Check column on both sides
    for i in range(n):
        if board[i][col] == 1:
            return False
    
    # Check diagonals in all directions
    for i in range(n):
        for j in range(n):
            if board[i][j] == 1:
                if abs(row - i) == abs(col - j):
                    return False
    
    return True

def solve_n_queens(board, queens_placed, n, forbidden):
    # Base case: If all queens are placed
    if queens_placed == n:
        return True
    
    # Try placing queen in each possible position
    for col in range(n):
        for row in range(n):
            # Skip if position is forbidden or already has a queen
            if (row, col) in forbidden or board[row][col] == 1:
                continue
                
            # Check if queen can be placed
            if is_safe(board, row, col, n):
                board[row][col] = 1
                
                # Recur to place rest of the queens
                if solve_n_queens(board, queens_placed + 1, n, forbidden):
                    return True
                
                # If placing queen doesn't lead to a solution, backtrack
                board[row][col] = 0
    
    return False

# Initialize the board
n = 8
board = [[0 for x in range(n)] for y in range(n)]
# Place the existing queen
board[0][7] = 1
# Mark forbidden position
forbidden = {(4,1)}

# Solve the puzzle starting with 1 queen already placed
if solve_n_queens(board, 1, n, forbidden):
    # Collect queen positions
    positions = []
    for i in range(n):
        for j in range(n):
            if board[i][j] == 1:
                positions.append(f"{i} {j}")
    print(", ".join(positions))
else:
    print("No solution exists")