def is_safe(board, row, col, n):
    # Check row
    for j in range(n):
        if board[row][j] == 1:
            return False
    
    # Check column
    for i in range(n):
        if board[i][col] == 1:
            return False
    
    # Check diagonal (top-left to bottom-right)
    for i, j in zip(range(row, -1, -1), range(col, -1, -1)):
        if board[i][j] == 1:
            return False
    for i, j in zip(range(row, n, 1), range(col, n, 1)):
        if board[i][j] == 1:
            return False
    
    # Check diagonal (top-right to bottom-left)
    for i, j in zip(range(row, -1, -1), range(col, n, 1)):
        if board[i][j] == 1:
            return False
    for i, j in zip(range(row, n, 1), range(col, -1, -1)):
        if board[i][j] == 1:
            return False
    
    return True

def solve_n_queens(board, n):
    def solve(row):
        if row >= n:
            return True
        
        # Skip rows that already have queens
        while row < n and 1 in board[row]:
            row += 1
        if row >= n:
            return True
        
        for col in range(n):
            # Skip forbidden positions
            if board[row][col] == 'X':
                continue
                
            if is_safe(board, row, col, n):
                board[row][col] = 1
                if solve(row + 1):
                    return True
                board[row][col] = 0
        return False
    
    if solve(0):
        positions = []
        for i in range(n):
            for j in range(n):
                if board[i][j] == 1:
                    positions.append(f"{i} {j}")
        return ", ".join(positions)
    return "No solution exists"

# Initialize the board
board = [
    [0,0,0,0,0,0,0,0],
    [1,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0],
    [0,0,0,0,'X',0,0,0],
    [0,0,0,0,0,0,0,0],
    [0,0,0,0,'X',0,0,'X'],
    [0,0,0,0,0,0,0,0],
    [0,0,0,1,0,0,0,0]
]

print(solve_n_queens(board, 8))