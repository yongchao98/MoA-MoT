def is_safe(board, row, col, n):
    # Check row
    for j in range(n):
        if board[row][j] == 1:
            return False
    
    # Check column
    for i in range(n):
        if board[i][col] == 1:
            return False
    
    # Check upper right diagonal
    for i, j in zip(range(row-1, -1, -1), range(col+1, n)):
        if board[i][j] == 1:
            return False
    
    # Check upper left diagonal
    for i, j in zip(range(row-1, -1, -1), range(col-1, -1, -1)):
        if board[i][j] == 1:
            return False
            
    # Check lower right diagonal
    for i, j in zip(range(row+1, n), range(col+1, n)):
        if board[i][j] == 1:
            return False
            
    # Check lower left diagonal
    for i, j in zip(range(row+1, n), range(col-1, -1, -1)):
        if board[i][j] == 1:
            return False
    
    return True

def solve_n_queens(board, row, n, forbidden, queens_placed):
    if queens_placed == n:
        return True
    
    if row >= n:
        return False
    
    for col in range(n):
        # Skip if position is forbidden or in the same row as fixed queen
        if (row, col) in forbidden:
            continue
            
        if is_safe(board, row, col, n):
            board[row][col] = 1
            
            if solve_n_queens(board, row + 1, n, forbidden, queens_placed + 1):
                return True
                
            board[row][col] = 0
    
    return False

def get_queen_positions(board):
    positions = []
    for i in range(len(board)):
        for j in range(len(board)):
            if board[i][j] == 1:
                positions.append(f"{i} {j}")
    return positions

# Initialize the board
n = 8
board = [[0 for x in range(n)] for y in range(n)]
board[7][5] = 1  # Given queen

# Forbidden positions (X)
forbidden = {(0,5), (1,7), (2,2), (6,3)}

# Start with one queen already placed
queens_placed = 1

# Solve the puzzle
if solve_n_queens(board, 0, n, forbidden, queens_placed):
    positions = get_queen_positions(board)
    print("<<<" + ", ".join(positions) + ">>>")
else:
    print("No solution exists")