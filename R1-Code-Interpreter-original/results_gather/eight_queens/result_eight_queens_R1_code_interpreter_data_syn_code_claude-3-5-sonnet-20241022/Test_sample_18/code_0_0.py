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

def solve_n_queens(n=8):
    # Initialize board with forbidden positions
    board = [[0 for x in range(n)] for y in range(n)]
    # Mark forbidden positions with -1
    forbidden = [(2,6), (3,5), (5,6), (7,0)]
    for row, col in forbidden:
        board[row][col] = -1
    
    def solve_util(board, col):
        if col >= n:
            return True
        
        for i in range(n):
            if board[i][col] == 0 and is_safe(board, i, col, n):
                board[i][col] = 1
                if solve_util(board, col + 1):
                    return True
                board[i][col] = 0
        return False
    
    if solve_util(board, 0) == False:
        print("Solution does not exist")
        return None
    
    # Get queen positions
    queens = []
    for i in range(n):
        for j in range(n):
            if board[i][j] == 1:
                queens.append(f"{i} {j}")
    
    return ", ".join(queens)

result = solve_n_queens()
print(result)