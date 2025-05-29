def create_board():
    board = [['.'] * 13 for _ in range(13)]
    
    # White stones (O)
    white = [
        (7,12), # H13
        (1,10), # B11
        (3,9),  # D10
        (3,6),  # D7
        (9,6),  # J7
        (10,5), # K6
        (11,5), # L6
        (8,4),  # I5
        (10,4), # K5
        (11,2)  # L3
    ]
    
    # Black stones (X)
    black = [
        (11,9), # L10
        (4,8),  # E9
        (4,7),  # E8
        (10,6), # K7
        (11,6), # L7
        (1,5),  # B6
        (9,5),  # J6
        (12,5), # M6
        (9,4),  # J5
        (10,3), # K4
        (5,1)   # F2
    ]
    
    for x, y in white:
        board[13-y][x] = 'O'
    for x, y in black:
        board[13-y][x] = 'X'
    
    return board

def check_group_liberties(board, x, y, visited=None):
    if visited is None:
        visited = set()
    
    if x < 0 or x >= 13 or y < 0 or y >= 13:
        return set()
    
    if board[y][x] == '.':
        return {(x,y)}
    
    if (x,y) in visited:
        return set()
    
    if board[y][x] not in ['O', 'X']:
        return set()
        
    visited.add((x,y))
    stone_type = board[y][x]
    liberties = set()
    
    for dx, dy in [(0,1), (1,0), (0,-1), (-1,0)]:
        nx, ny = x + dx, y + dy
        if 0 <= nx < 13 and 0 <= ny < 13:
            if board[ny][nx] == '.':
                liberties.add((nx,ny))
            elif board[ny][nx] == stone_type:
                liberties.update(check_group_liberties(board, nx, ny, visited))
    
    return liberties

def would_capture(board, x, y):
    if board[y][x] != '.':
        return 0
    
    # Try placing black stone
    board[y][x] = 'X'
    captured = 0
    
    # Check adjacent white groups
    for dx, dy in [(0,1), (1,0), (0,-1), (-1,0)]:
        nx, ny = x + dx, y + dy
        if 0 <= nx < 13 and 0 <= ny < 13 and board[ny][nx] == 'O':
            if not check_group_liberties(board, nx, ny):
                captured += 1
    
    # Restore board
    board[y][x] = '.'
    return captured

def find_best_capture():
    board = create_board()
    max_captures = 0
    best_move = None
    
    for y in range(13):
        for x in range(13):
            if board[y][x] == '.':
                captures = would_capture(board, x, y)
                if captures > max_captures:
                    max_captures = captures
                    # Convert to Go coordinates
                    col = chr(x + ord('A'))
                    row = 13 - y
                    best_move = (col, row)
    
    if best_move:
        print(f"<<<{best_move[0]}{best_move[1]}>>>")
    else:
        print("<<<No capturing move found>>>")

find_best_capture()