def simulate_move(board, move):
    # Convert move to board indices
    col, row = move
    col_idx = ord(col) - ord('A')
    row_idx = 9 - int(row)
    
    # Place the black stone
    board[row_idx][col_idx] = 'X'
    
    # Check for captures
    def capture_group(r, c, color, visited):
        if r < 0 or r >= 9 or c < 0 or c >= 9:
            return False
        if board[r][c] == '.' or board[r][c] == 'X':
            return True
        if (r, c) in visited:
            return True
        if board[r][c] != color:
            return False
        
        visited.add((r, c))
        up = capture_group(r-1, c, color, visited)
        down = capture_group(r+1, c, color, visited)
        left = capture_group(r, c-1, color, visited)
        right = capture_group(r, c+1, color, visited)
        
        return up and down and left and right
    
    captured_stones = 0
    for r in range(9):
        for c in range(9):
            if board[r][c] == 'O':
                visited = set()
                if capture_group(r, c, 'O', visited):
                    captured_stones += len(visited)
                    for vr, vc in visited:
                        board[vr][vc] = '.'
    
    return captured_stones

# Initial board setup
initial_board = [
    ['.', 'X', '.', '.', '.', '.', '.', '.', '.'],
    ['O', '.', 'X', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', 'X', '.', '.', '.', '.', '.'],
    ['.', '.', 'X', '.', 'X', '.', '.', 'X', '.'],
    ['.', 'X', 'O', 'O', 'O', 'X', '.', '.', '.'],
    ['O', '.', 'X', 'O', 'X', '.', '.', '.', '.'],
    ['O', '.', '.', 'X', '.', '.', '.', '.', '.'],
    ['.', '.', '.', 'O', '.', 'O', '.', '.', '.'],
    ['.', '.', 'O', '.', '.', 'X', '.', '.', '.']
]

# Test potential moves
moves = ['C5', 'E4']
results = {move: simulate_move([row[:] for row in initial_board], move) for move in moves}

# Output the results
print(results)