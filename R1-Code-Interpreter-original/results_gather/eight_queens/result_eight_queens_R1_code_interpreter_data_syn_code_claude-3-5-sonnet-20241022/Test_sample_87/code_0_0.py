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
    
    for i in range(n):
        if (i, col) not in forbidden and is_safe(board, i, col, n):
            board[i][col] = 1
            if solve_n_queens(board, col + 1, n, forbidden):
                return True
            board[i][col] = 0
    
    return False

def get_queen_positions(board):
    positions = []
    for i in range(len(board)):
        for j in range(len(board)):
            if board[i][j] == 1:
                positions.append(f"{i} {j}")
    return positions

# Initialize the board and forbidden positions
n = 8
board = [[0 for x in range(n)] for y in range(n)]
forbidden = {(1, 4), (2, 6), (5, 5)}  # X positions

# Solve the puzzle
if solve_n_queens(board, 0, n, forbidden):
    positions = get_queen_positions(board)
    print("<<<" + ", ".join(positions) + ">>>")
else:
    print("No solution exists")