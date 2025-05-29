def is_safe(board, row, col, cols, diag1, diag2):
    if cols[col] or diag1[row - col] or diag2[row + col] or board[row][col] == 'X':
        return False
    return True

def solve_n_queens(board, row, cols, diag1, diag2, positions):
    if row == 8:
        return True
    
    for col in range(8):
        if is_safe(board, row, col, cols, diag1, diag2):
            # Place the queen
            cols[col] = True
            diag1[row - col] = True
            diag2[row + col] = True
            positions.append((row, col))
            
            # Recur to place the rest
            if solve_n_queens(board, row + 1, cols, diag1, diag2, positions):
                return True
            
            # Backtrack
            cols[col] = False
            diag1[row - col] = False
            diag2[row + col] = False
            positions.pop()
    
    return False

def find_queen_positions():
    board = [
        ['0', '0', '0', '0', '0', '0', '0', '0'],
        ['0', '0', '0', '0', '0', '0', '0', '0'],
        ['0', '0', '0', '0', '0', '0', '0', '0'],
        ['0', '0', '0', '0', '0', '0', 'X', '0'],
        ['0', '0', '0', '0', '0', '0', '0', '0'],
        ['0', '0', '0', '0', '0', '0', '0', '0'],
        ['0', '0', '0', '0', '0', '0', '0', '0'],
        ['0', '0', '0', '0', '0', '0', '0', '0']
    ]
    
    cols = [False] * 8
    diag1 = [False] * 15
    diag2 = [False] * 15
    positions = []
    
    if solve_n_queens(board, 0, cols, diag1, diag2, positions):
        return positions
    else:
        return []

positions = find_queen_positions()
print(positions)