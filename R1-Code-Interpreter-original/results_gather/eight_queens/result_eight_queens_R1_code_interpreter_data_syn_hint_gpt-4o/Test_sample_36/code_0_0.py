def is_safe(board, row, col, n, cols, diag1, diag2):
    if cols[col] or diag1[row - col] or diag2[row + col] or board[row][col] == 'X':
        return False
    return True

def solve_n_queens(board, row, n, cols, diag1, diag2, positions):
    if row == n:
        return True
    
    for col in range(n):
        if is_safe(board, row, col, n, cols, diag1, diag2):
            # Place the queen
            cols[col] = True
            diag1[row - col] = True
            diag2[row + col] = True
            positions.append((row, col))
            
            # Recur to place the rest
            if solve_n_queens(board, row + 1, n, cols, diag1, diag2, positions):
                return True
            
            # Backtrack
            cols[col] = False
            diag1[row - col] = False
            diag2[row + col] = False
            positions.pop()
    
    return False

def find_queen_positions():
    n = 8
    board = [
        ['0', '0', '0', '0', '0', '0', '0', '0'],
        ['0', '0', '0', '0', '0', '0', '0', '0'],
        ['0', '0', '0', 'X', '0', '0', '0', '0'],
        ['0', '0', '0', '0', '0', '0', '0', 'X'],
        ['0', '0', '0', '0', '0', '0', '0', '0'],
        ['0', '0', '0', '0', '0', '0', '0', 'X'],
        ['0', '0', '0', '0', '0', '0', '0', '0'],
        ['0', '0', '0', '0', 'X', '0', '0', '0']
    ]
    
    cols = [False] * n
    diag1 = [False] * (2 * n)
    diag2 = [False] * (2 * n)
    positions = []
    
    if solve_n_queens(board, 0, n, cols, diag1, diag2, positions):
        return positions
    else:
        return []

queen_positions = find_queen_positions()
print(queen_positions)