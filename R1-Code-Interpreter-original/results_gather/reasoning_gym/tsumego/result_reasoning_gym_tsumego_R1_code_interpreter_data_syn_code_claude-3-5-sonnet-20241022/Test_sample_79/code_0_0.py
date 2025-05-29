def is_valid_pos(x, y):
    return 0 <= x < 10 and 0 <= y < 10

def get_neighbors(x, y):
    return [(x+1, y), (x-1, y), (x, y+1), (x, y-1)]

def count_captures(board, move_x, move_y):
    if board[move_y][move_x] != '.':
        return 0
        
    # Temporarily place the black stone
    board[move_y][move_x] = 'X'
    
    captured = set()
    checked = set()
    
    def check_group(x, y, stone_type):
        if (x, y) in checked:
            return False
        if not is_valid_pos(x, y):
            return False
        if board[y][x] == '.':
            return False
        if board[y][x] != stone_type:
            return True
            
        group = set([(x, y)])
        liberty = False
        to_check = [(x, y)]
        
        while to_check:
            curr_x, curr_y = to_check.pop()
            checked.add((curr_x, curr_y))
            
            for nx, ny in get_neighbors(curr_x, curr_y):
                if not is_valid_pos(nx, ny):
                    continue
                if board[ny][nx] == '.':
                    liberty = True
                elif board[ny][nx] == stone_type and (nx, ny) not in group:
                    group.add((nx, ny))
                    to_check.append((nx, ny))
                    
        if not liberty:
            captured.update(group)
        return liberty

    # Check all adjacent positions for white groups
    for nx, ny in get_neighbors(move_x, move_y):
        if is_valid_pos(nx, ny) and board[ny][nx] == 'O':
            check_group(nx, ny, 'O')
    
    # Restore the board
    board[move_y][move_x] = '.'
    
    return len(captured)

# Initialize the board
board = [
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', 'X', 'O', '.'],
    ['O', '.', '.', '.', '.', '.', '.', 'X', 'O', 'O'],
    ['.', '.', '.', '.', '.', '.', 'X', 'O', 'X', 'X'],
    ['.', '.', '.', '.', 'X', '.', '.', '.', 'X', 'O'],
    ['.', '.', '.', 'O', '.', '.', '.', 'X', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', 'X', 'O', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', 'X', '.', '.'],
    ['.', '.', '.', '.', '.', 'X', '.', '.', '.', '.']
]

# Try all possible moves and find the one with maximum captures
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
    # Convert to Go coordinates (A-J, 1-10)
    col = chr(ord('A') + best_move[0])
    row = 10 - best_move[1]
    print(f"Best move: {col}{row} capturing {max_captures} stones")