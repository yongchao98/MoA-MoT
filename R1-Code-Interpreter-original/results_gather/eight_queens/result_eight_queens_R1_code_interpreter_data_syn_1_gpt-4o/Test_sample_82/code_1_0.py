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

def find_queen_positions():
    board = [
        ['0', '0', '1', '0', '0', '0', '0', '0'],
        ['0', '0', '0', '0', '0', '0', '0', '0'],
        ['0', '0', '0', '0', '0', '0', '0', '0'],
        ['0', 'X', '0', '0', '0', '0', '0', '0'],
        ['0', '0', '0', 'X', '0', '0', '0', '0'],
        ['0', '0', '0', '0', '0', '0', '0', '0'],
        ['0', '0', '0', '0', '0', '0', '0', '0'],
        ['0', '0', '0', '0', '0', '0', 'X', '0']
    ]
    
    n = 8
    cols = [False] * n
    diag1 = [False] * (2 * n)
    diag2 = [False] * (2 * n)
    positions = []
    
    # Mark the initial queens
    for r in range(n):
        for c in range(n):
            if board[r][c] == '1':
                cols[c] = True
                diag1[r - c] = True
                diag2[r + c] = True
                positions.append((r, c))
    
    # Start from the first row without a queen
    start_row = 0
    for r in range(n):
        if any(board[r][c] == '1' for c in range(n)):
            start_row = r + 1
    
    solve_n_queens(board, start_row, n, cols, diag1, diag2, positions)
    
    return positions

queen_positions = find_queen_positions()
print(queen_positions)