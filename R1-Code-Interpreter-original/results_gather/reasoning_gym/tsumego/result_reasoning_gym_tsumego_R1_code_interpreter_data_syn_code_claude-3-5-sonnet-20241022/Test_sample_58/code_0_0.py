def get_neighbors(x, y):
    return [(x+1,y), (x-1,y), (x,y+1), (x,y-1)]

def count_captures(board, move_x, move_y):
    if board[move_y][move_x] != '.':
        return 0
        
    # Temporarily place black stone
    board[move_y][move_x] = 'X'
    captured = set()
    checked = set()
    
    # Check each neighbor for white groups with no liberties
    for nx, ny in get_neighbors(move_x, move_y):
        if 0 <= nx < 10 and 0 <= ny < 10 and board[ny][nx] == 'O' and (nx,ny) not in checked:
            group = set()
            liberties = set()
            stack = [(nx,ny)]
            
            while stack:
                cx, cy = stack.pop()
                if (cx,cy) in group:
                    continue
                group.add((cx,cy))
                
                for adj_x, adj_y in get_neighbors(cx, cy):
                    if 0 <= adj_x < 10 and 0 <= adj_y < 10:
                        if board[adj_y][adj_x] == 'O':
                            stack.append((adj_x, adj_y))
                        elif board[adj_y][adj_x] == '.':
                            liberties.add((adj_x, adj_y))
            
            checked.update(group)
            if not liberties:
                captured.update(group)
    
    # Restore the board
    board[move_y][move_x] = '.'
    return len(captured)

# Initialize board
board = [
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', 'X', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', 'X', '.', '.', '.', '.'],
    ['.', 'O', '.', '.', '.', '.', '.', 'O', '.', '.'],
    ['X', 'O', 'X', 'X', 'X', '.', '.', '.', '.', '.'],
    ['O', 'O', 'O', 'X', '.', '.', '.', '.', '.', '.'],
    ['O', '.', 'X', '.', '.', '.', '.', '.', '.', '.'],
    ['X', 'O', 'O', '.', '.', '.', '.', 'O', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.']
]

# Find best move
max_captures = 0
best_move = None

for y in range(10):
    for x in range(10):
        captures = count_captures(board, x, y)
        if captures > max_captures:
            max_captures = captures
            best_move = (x, y)

if best_move:
    col = chr(ord('A') + best_move[0])
    row = 10 - best_move[1]
    print(f"Best move is {col}{row} capturing {max_captures} stones")