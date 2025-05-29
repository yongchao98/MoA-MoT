def is_safe(board, row, col, n):
    # Check if position is forbidden (X)
    if row == 4 and col == 1:
        return False
    
    # Check row on left side
    for j in range(col):
        if board[row][j] == 1:
            return False
    
    # Check upper diagonal on left side
    for i, j in zip(range(row, -1, -1), range(col, -1, -1)):
        if board[i][j] == 1:
            return False
    
    # Check lower diagonal on left side
    for i, j in zip(range(row, n, 1), range(col, -1, -1)):
        if board[i][j] == 1:
            return False
    
    return True

def solve_n_queens(board, col, n):
    if col >= n:
        return True
    
    for i in range(n):
        if is_safe(board, i, col, n):
            board[i][col] = 1
            
            if solve_n_queens(board, col + 1, n):
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

# Initialize board
n = 8
board = [[0 for x in range(n)] for y in range(n)]

# Solve the puzzle
if solve_n_queens(board, 0, n):
    print_solution(board)
else:
    print("No solution exists")