def is_valid(x, y):
    return 0 <= x < 9 and 0 <= y < 9

def get_neighbors(x, y):
    return [(x+dx, y+dy) for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)] if is_valid(x+dx, y+dy)]

def count_captures(board, x, y, stone):
    if board[y][x] != '.':
        return 0
    
    # Place stone temporarily
    board[y][x] = stone
    captured = set()
    checked = set()
    
    # Check all neighbors for opponent groups that can be captured
    for nx, ny in get_neighbors(x, y):
        if (nx, ny) not in checked and board[ny][nx] == 'O':  # opponent is white
            group = set()
            liberties = set()
            stack = [(nx, ny)]
            
            while stack:
                cx, cy = stack.pop()
                if (cx, cy) not in group:
                    group.add((cx, cy))
                    for ax, ay in get_neighbors(cx, cy):
                        if board[ay][ax] == 'O':
                            stack.append((ax, ay))
                        elif board[ay][ax] == '.':
                            liberties.add((ax, ay))
            
            checked.update(group)
            if not liberties:
                captured.update(group)
    
    # Remove temporary stone
    board[y][x] = '.'
    
    return len(captured)

# Initialize board
board = [
    ['.', '.', '.', '.', 'X', '.', '.', 'X', 'X'],
    ['.', '.', '.', '.', '.', 'X', 'X', 'O', 'O'],
    ['.', '.', '.', '.', '.', '.', 'X', 'O', '.'],
    ['.', 'O', 'O', '.', '.', '.', 'X', 'X', 'O'],
    ['.', 'O', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', 'X', '.', '.', '.', '.', 'O'],
    ['.', '.', '.', 'X', '.', '.', '.', '.', 'X'],
    ['.', '.', '.', '.', '.', '.', 'X', '.', '.']
]

# Try all empty positions
max_captures = 0
best_move = None

for y in range(9):
    for x in range(9):
        if board[y][x] == '.':
            captures = count_captures(board, x, y, 'X')
            if captures > max_captures:
                max_captures = captures
                best_move = (x, y)

if best_move:
    col = chr(ord('A') + best_move[0])
    row = 9 - best_move[1]
    print(f"Best move is {col}{row} capturing {max_captures} stones")