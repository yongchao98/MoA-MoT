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

def solve_n_queens(board, col, n, forbidden):
    if col >= n:
        return True
    
    # Check if this column already has a queen
    queen_in_col = False
    for i in range(n):
        if board[i][col] == 1:
            queen_in_col = True
            break
    
    if queen_in_col:
        return solve_n_queens(board, col + 1, n, forbidden)
    
    for i in range(n):
        if (i, col) not in forbidden and is_safe(board, i, col, n):
            board[i][col] = 1
            if solve_n_queens(board, col + 1, n, forbidden):
                return True
            board[i][col] = 0
    
    return False

# Initialize the board and forbidden positions
n = 8
board = [[0 for x in range(n)] for y in range(n)]
forbidden = {(1, 4), (6, 2)}  # X positions

# Try to solve
if solve_n_queens(board, 0, n, forbidden):
    # Convert solution to required format
    queens = []
    for i in range(n):
        for j in range(n):
            if board[i][j] == 1:
                queens.append(f"{i} {j}")
    print("<<<" + ", ".join(queens) + ">>>")
else:
    print("No solution exists")