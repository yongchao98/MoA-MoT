def is_safe(board, row, col, n=8):
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

def solve_n_queens(forbidden_positions):
    n = 8
    board = [[0 for x in range(n)] for y in range(n)]
    
    # Mark forbidden positions with -1
    for pos in forbidden_positions:
        board[pos[0]][pos[1]] = -1
    
    def solve_util(board, col):
        if col >= n:
            return True
        
        for i in range(n):
            if board[i][col] == 0 and is_safe(board, i, col):
                board[i][col] = 1
                
                if solve_util(board, col + 1):
                    return True
                
                board[i][col] = 0
        
        return False
    
    if solve_util(board, 0) == False:
        print("Solution does not exist")
        return None
    
    # Collect queen positions
    queen_positions = []
    for i in range(n):
        for j in range(n):
            if board[i][j] == 1:
                queen_positions.append(f"{i} {j}")
    
    return ", ".join(queen_positions)

# Define forbidden positions
forbidden = [(1,3), (2,7), (5,5)]

result = solve_n_queens(forbidden)
print(result)