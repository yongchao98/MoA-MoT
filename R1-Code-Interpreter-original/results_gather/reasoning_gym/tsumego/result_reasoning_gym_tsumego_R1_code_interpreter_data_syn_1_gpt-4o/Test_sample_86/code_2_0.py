# Simulate the move at I9 and check the capture
def simulate_move(board, move):
    # Convert move to board indices
    col, row = move
    col_index = ord(col) - ord('A')
    row_index = 10 - int(row)
    
    # Place the black stone
    board[row_index][col_index] = 'X'
    
    # Check for captures
    captures = []
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    for dr, dc in directions:
        r, c = row_index + dr, col_index + dc
        if 0 <= r < 10 and 0 <= c < 10 and board[r][c] == 'O':
            # Check if this white stone is captured
            if is_captured(board, r, c):
                captures.append((r, c))
    
    return captures

def is_captured(board, r, c):
    # Check if the stone at (r, c) is captured
    visited = set()
    return not has_liberty(board, r, c, visited)

def has_liberty(board, r, c, visited):
    if (r, c) in visited:
        return False
    visited.add((r, c))
    
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    for dr, dc in directions:
        nr, nc = r + dr, c + dc
        if 0 <= nr < 10 and 0 <= nc < 10:
            if board[nr][nc] == '.':
                return True
            if board[nr][nc] == 'O' and has_liberty(board, nr, nc, visited):
                return True
    return False

# Initial board setup
board = [
    ['.', '.', '.', '.', '.', '.', '.', '.', 'O', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'O'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['X', 'O', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', 'O', '.', '.', '.', '.'],
    ['.', '.', '.', '.', 'O', '.', '.', 'X', '.', '.'],
    ['.', '.', '.', '.', '.', '.', 'X', 'O', 'X', '.'],
    ['O', '.', 'O', '.', '.', 'X', 'O', 'O', '.', 'X'],
    ['.', '.', 'O', '.', '.', 'O', 'X', 'O', 'X', '.']
]

# Simulate the move at I9
captures = simulate_move(board, ('I', '9'))
print(captures)