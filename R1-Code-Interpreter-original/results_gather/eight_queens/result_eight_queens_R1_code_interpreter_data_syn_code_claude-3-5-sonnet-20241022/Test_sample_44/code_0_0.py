def is_safe(board, row, col, n):
    # Check row
    for j in range(n):
        if board[row][j] == 1:
            return False
    
    # Check column
    for i in range(n):
        if board[i][col] == 1:
            return False
    
    # Check upper left diagonal
    for i, j in zip(range(row, -1, -1), range(col, -1, -1)):
        if board[i][j] == 1:
            return False
    
    # Check upper right diagonal
    for i, j in zip(range(row, -1, -1), range(col, n)):
        if board[i][j] == 1:
            return False
    
    # Check lower left diagonal
    for i, j in zip(range(row, n), range(col, -1, -1)):
        if board[i][j] == 1:
            return False
    
    # Check lower right diagonal
    for i, j in zip(range(row, n), range(col, n)):
        if board[i][j] == 1:
            return False
    
    return True

def solve_n_queens(board, n):
    # First, count existing queens and store their positions
    queen_positions = []
    queen_count = 0
    for i in range(n):
        for j in range(n):
            if board[i][j] == 1:
                queen_count += 1
                queen_positions.append((i, j))
    
    def solve_util(row):
        nonlocal queen_count, queen_positions
        
        if queen_count == 8:
            return True
        
        if row >= n:
            return False
        
        for col in range(n):
            if board[row][col] != 'X' and board[row][col] != 1:
                if is_safe(board, row, col, n):
                    board[row][col] = 1
                    queen_count += 1
                    queen_positions.append((row, col))
                    
                    if solve_util(row + 1):
                        return True
                    
                    board[row][col] = 0
                    queen_count -= 1
                    queen_positions.pop()
        
        return solve_util(row + 1)
    
    solve_util(0)
    return queen_positions

# Initialize the board
board = [
    [0,0,0,0,0,'X',0,0],
    [0,0,0,0,0,0,0,'X'],
    [0,'X',0,0,1,0,0,0],
    [0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0],
    [0,0,0,0,0,1,0,0],
    [0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0]
]

positions = solve_n_queens(board, 8)
positions.sort()  # Sort positions for consistent output
result = ", ".join(f"{pos[0]} {pos[1]}" for pos in positions)
print(f"<<<{result}>>>")