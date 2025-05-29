def is_safe(board, row, col, n):
    # Check row on left side
    for j in range(col):
        if board[row][j] == 1:
            return False
    
    # Check upper diagonal on left side
    for i, j in zip(range(row-1, -1, -1), range(col-1, -1, -1)):
        if board[i][j] == 1:
            return False
    
    # Check lower diagonal on left side
    for i, j in zip(range(row+1, n), range(col-1, -1, -1)):
        if board[i][j] == 1:
            return False
    
    return True

def solve_queens(board, col, n, forbidden):
    if col >= n:
        return True
    
    # Skip if we already have a queen in this column
    queen_in_col = False
    for i in range(n):
        if board[i][col] == 1:
            queen_in_col = True
            break
    
    if queen_in_col:
        return solve_queens(board, col + 1, n, forbidden)
    
    for i in range(n):
        if (i, col) == forbidden:  # Skip forbidden position
            continue
            
        if is_safe(board, i, col, n):
            board[i][col] = 1
            
            if solve_queens(board, col + 1, n, forbidden):
                return True
                
            board[i][col] = 0
    
    return False

def print_solution(board):
    positions = []
    for i in range(len(board)):
        for j in range(len(board)):
            if board[i][j] == 1:
                positions.append(f"{i} {j}")
    print("<<<" + ", ".join(positions) + ">>>")

# Initialize board with given queen and constraints
n = 8
board = [[0 for x in range(n)] for y in range(n)]
board[4][0] = 1  # Given queen
forbidden = (7, 5)  # Forbidden position

# Solve the puzzle
if solve_queens(board, 0, n, forbidden):
    print_solution(board)
else:
    print("No solution exists")