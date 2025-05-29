def is_safe(board, row, col, n, cols, diag1, diag2):
    if cols[col] or diag1[row - col] or diag2[row + col] or board[row][col] == 'X':
        return False
    return True

def solve_n_queens(board, row, n, cols, diag1, diag2, positions):
    if row == n:
        return True
    
    for col in range(n):
        if is_safe(board, row, col, n, cols, diag1, diag2):
            board[row][col] = '1'
            cols[col] = True
            diag1[row - col] = True
            diag2[row + col] = True
            positions.append((row, col))
            
            if solve_n_queens(board, row + 1, n, cols, diag1, diag2, positions):
                return True
            
            board[row][col] = '0'
            cols[col] = False
            diag1[row - col] = False
            diag2[row + col] = False
            positions.pop()
    
    return False

def place_queens():
    n = 8
    board = [['0'] * n for _ in range(n)]
    board[5][1] = 'X'
    
    cols = [False] * n
    diag1 = [False] * (2 * n)
    diag2 = [False] * (2 * n)
    
    positions = []
    
    solve_n_queens(board, 0, n, cols, diag1, diag2, positions)
    
    return positions

positions = place_queens()
print(positions)