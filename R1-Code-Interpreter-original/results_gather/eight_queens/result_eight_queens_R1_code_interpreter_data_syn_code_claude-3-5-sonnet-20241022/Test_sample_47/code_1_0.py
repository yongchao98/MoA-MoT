def is_safe(board, row, col, n):
    # Check row on both sides
    for j in range(n):
        if board[row][j] == 1:
            return False
    
    # Check column
    for i in range(n):
        if board[i][col] == 1:
            return False
    
    # Check diagonals
    for i in range(n):
        for j in range(n):
            if board[i][j] == 1:
                if abs(row - i) == abs(col - j):
                    return False
    
    return True

def solve_n_queens(board, queens_placed, n):
    if queens_placed == 8:
        return True

    for col in range(n):
        for row in range(n):
            # Skip if cell is blocked or already has a queen
            if board[row][col] == 2 or board[row][col] == 1:
                continue
                
            if is_safe(board, row, col, n):
                board[row][col] = 1
                
                if solve_n_queens(board, queens_placed + 1, n):
                    return True
                    
                board[row][col] = 0
    
    return False

# Initialize the board
n = 8
board = [
    [0,0,0,0,0,0,0,2],  # X marked as 2
    [0,0,0,1,0,0,0,0],
    [0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0],
    [2,0,0,2,0,0,0,0],
    [0,0,2,0,0,0,0,0],
    [0,0,0,0,0,0,0,0]
]

# Count initial queens
initial_queens = sum(row.count(1) for row in board)

# Solve the puzzle
if solve_n_queens(board, initial_queens, n):
    # Print positions of queens
    positions = []
    for i in range(n):
        for j in range(n):
            if board[i][j] == 1:
                positions.append(f"{i} {j}")
    print(", ".join(positions))
else:
    print("No solution exists")