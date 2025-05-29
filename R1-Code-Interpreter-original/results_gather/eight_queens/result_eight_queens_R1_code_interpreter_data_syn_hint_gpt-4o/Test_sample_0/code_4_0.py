def is_valid(board, row, col, cols, diag1, diag2):
    if col in cols or (row - col) in diag1 or (row + col) in diag2:
        return False
    return True

def solve_n_queens(board, row, positions, cols, diag1, diag2):
    if row == len(board):
        return positions
    
    for col in range(len(board)):
        if board[row][col] == 0 and is_valid(board, row, col, cols, diag1, diag2):
            board[row][col] = 1
            positions.append((row, col))
            cols.add(col)
            diag1.add(row - col)
            diag2.add(row + col)
            
            result = solve_n_queens(board, row + 1, positions, cols, diag1, diag2)
            if result:
                return result
            
            board[row][col] = 0
            positions.pop()
            cols.remove(col)
            diag1.remove(row - col)
            diag2.remove(row + col)
    
    return None

# Initial board setup
board = [
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, -1, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, -1],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0]
]

# Solve the problem
positions = solve_n_queens(board, 0, [], set(), set(), set())
print(positions)