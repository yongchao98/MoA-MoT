def is_valid_pos(x, y):
    return 0 <= x < 10 and 0 <= y < 10

def get_neighbors(x, y):
    return [(x+dx, y+dy) for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)] if is_valid_pos(x+dx, y+dy)]

def count_captures(board, move_x, move_y):
    if board[move_y][move_x] != '.':
        return 0
        
    # Temporarily place the black stone
    board[move_y][move_x] = 'X'
    captured = set()
    
    def check_group(x, y, color, visited):
        if (x, y) in visited:
            return True
        if not is_valid_pos(x, y) or board[y][x] == '.':
            return False
        if board[y][x] != color:
            return True
            
        visited.add((x, y))
        has_liberty = False
        
        for nx, ny in get_neighbors(x, y):
            if not check_group(nx, ny, color, visited):
                has_liberty = True
                
        return has_liberty

    # Check each neighboring white group
    for nx, ny in get_neighbors(move_x, move_y):
        if board[ny][nx] == 'O':
            group = set()
            if not check_group(nx, ny, 'O', group):
                captured.update(group)
    
    # Restore the board
    board[move_y][move_x] = '.'
    return len(captured)

# Initialize the board
board = [
    ['.', '.', '.', '.', 'O', '.', '.', 'X', 'O', 'X'],
    ['X', '.', 'X', '.', '.', '.', 'X', 'O', 'O', '.'],
    ['.', '.', '.', '.', '.', 'X', '.', 'O', 'O', 'X'],
    ['.', '.', '.', '.', '.', '.', '.', '.', 'X', 'X'],
    ['.', '.', 'X', '.', '.', '.', '.', '.', '.', 'O'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', 'O', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', 'X', '.', '.', '.'],
    ['.', '.', 'O', '.', '.', '.', '.', '.', '.', '.']
]

# Find the move with maximum captures
max_captures = 0
best_move = None

for y in range(10):
    for x in range(10):
        if board[y][x] == '.':
            captures = count_captures(board, x, y)
            if captures > max_captures:
                max_captures = captures
                best_move = (x, y)

if best_move:
    col = chr(ord('A') + best_move[0])
    row = 10 - best_move[1]
    print(f"Best move: {col}{row} with {max_captures} captures")