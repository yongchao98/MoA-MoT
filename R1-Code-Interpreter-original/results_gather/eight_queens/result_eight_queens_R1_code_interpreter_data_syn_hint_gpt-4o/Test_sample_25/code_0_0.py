def is_safe(board, row, col, cols, diag1, diag2):
    if col in cols or (row - col) in diag1 or (row + col) in diag2:
        return False
    return True

def solve_n_queens(board, row, cols, diag1, diag2, positions):
    if row == 8:
        return positions
    
    for col in range(8):
        if board[row][col] == 'X':
            continue
        if is_safe(board, row, col, cols, diag1, diag2):
            cols.add(col)
            diag1.add(row - col)
            diag2.add(row + col)
            positions.append((row, col))
            
            result = solve_n_queens(board, row + 1, cols, diag1, diag2, positions)
            if result:
                return result
            
            # Backtrack
            cols.remove(col)
            diag1.remove(row - col)
            diag2.remove(row + col)
            positions.pop()
    
    return None

# Initial board setup
board = [
    ['0', '0', '0', '0', '0', '0', '0', '0'],
    ['0', '0', '0', '0', '0', '0', '0', '0'],
    ['0', '0', 'X', '0', '0', '0', '0', '0'],
    ['0', 'X', '0', '0', 'X', '0', '0', '0'],
    ['0', '0', '0', '0', '0', '0', '0', '0'],
    ['0', '0', '0', '0', '0', '0', '0', '0'],
    ['0', '0', '0', '0', '0', '0', '0', '0'],
    ['0', '0', '0', '0', '0', '0', '0', 'X']
]

# Solve the problem
positions = solve_n_queens(board, 0, set(), set(), set(), [])
print(positions)