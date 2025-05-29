def create_board():
    # Create empty board
    board = [['.'] * 12 for _ in range(12)]
    
    # Place stones exactly as shown in the diagram
    # Black stones (X)
    black_positions = {
        'B12': (0,1), 'K10': (2,10), 'J10': (2,9),
        'L9': (3,11), 'I9': (3,8), 'I8': (4,8),
        'J7': (5,9), 'H6': (6,7), 'D6': (6,3),
        'K4': (8,10), 'C2': (10,2)
    }
    
    # White stones (O)
    white_positions = {
        'A11': (1,0), 'K11': (1,10), 'I10': (2,8),
        'J9': (3,9), 'K9': (3,10), 'O8': (4,3),
        'L8': (4,11), 'A3': (9,0), 'D2': (10,3)
    }
    
    # Place the stones
    for _, (x, y) in black_positions.items():
        board[x][y] = 'X'
    for _, (x, y) in white_positions.items():
        board[x][y] = 'O'
    
    return board

def get_group(board, x, y, color, visited=None):
    if visited is None:
        visited = set()
    
    if (x < 0 or x >= 12 or y < 0 or y >= 12 or 
        board[x][y] != color or (x,y) in visited):
        return set()
    
    group = {(x,y)}
    visited.add((x,y))
    
    for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
        new_x, new_y = x + dx, y + dy
        group.update(get_group(board, new_x, new_y, color, visited))
    
    return group

def get_liberties(board, group):
    liberties = set()
    for x, y in group:
        for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
            new_x, new_y = x + dx, y + dy
            if (0 <= new_x < 12 and 0 <= new_y < 12 and 
                board[new_x][new_y] == '.'):
                liberties.add((new_x, new_y))
    return liberties

def evaluate_move(board, x, y):
    if board[x][y] != '.':
        return 0
    
    # Make move
    board[x][y] = 'X'
    captured = 0
    
    # Check all adjacent positions for white groups
    checked = set()
    for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
        new_x, new_y = x + dx, y + dy
        if (0 <= new_x < 12 and 0 <= new_y < 12 and 
            board[new_x][new_y] == 'O' and 
            (new_x, new_y) not in checked):
            group = get_group(board, new_x, new_y, 'O')
            checked.update(group)
            liberties = get_liberties(board, group)
            if not liberties:
                captured += len(group)
    
    # Undo move
    board[x][y] = '.'
    return captured

def find_best_capture():
    board = create_board()
    best_move = None
    max_captures = 0
    
    for x in range(12):
        for y in range(12):
            captures = evaluate_move(board, x, y)
            if captures > max_captures:
                max_captures = captures
                # Convert to board coordinates
                col = chr(ord('A') + y)
                row = 12 - x
                best_move = f"{col}{row}"
    
    return best_move, max_captures

move, captures = find_best_capture()
print(f"Best move: {move} (captures {captures} stones)")