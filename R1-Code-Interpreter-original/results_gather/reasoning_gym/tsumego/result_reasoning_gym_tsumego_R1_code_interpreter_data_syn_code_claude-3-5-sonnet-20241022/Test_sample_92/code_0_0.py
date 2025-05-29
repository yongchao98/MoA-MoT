def create_board():
    board = [['.'] * 11 for _ in range(11)]
    # Add black stones (X)
    black_positions = [
        (6,0), (5,0), (3,0),  # left side
        (6,2),  # C5
        (0,6), (1,6), (2,6),  # G column
        (0,8),  # I11
        (3,8), (4,8),  # I positions
        (1,9), (0,10)  # bottom right
    ]
    # Add white stones (O)
    white_positions = [
        (5,4),  # E5
        (2,8), (2,9),  # I2, I3
        (1,8), (1,9),  # H2, H3
        (1,5)   # F1
    ]
    
    for x, y in black_positions:
        board[x][y] = 'X'
    for x, y in white_positions:
        board[x][y] = 'O'
    
    return board

def is_valid(x, y):
    return 0 <= x < 11 and 0 <= y < 11

def get_neighbors(x, y):
    directions = [(0,1), (1,0), (0,-1), (-1,0)]
    return [(x+dx, y+dy) for dx, dy in directions if is_valid(x+dx, y+dy)]

def count_captures(board, x, y):
    if board[x][y] != '.':
        return 0
    
    # Temporarily place stone
    board[x][y] = 'X'
    captured = set()
    
    # Check each direction for captures
    for nx, ny in get_neighbors(x, y):
        if board[nx][ny] == 'O':
            group = [(nx, ny)]
            visited = {(nx, ny)}
            liberties = set()
            
            # Find the connected group and its liberties
            while group:
                cx, cy = group.pop()
                for ax, ay in get_neighbors(cx, cy):
                    if (ax, ay) not in visited:
                        if board[ax][ay] == 'O':
                            group.append((ax, ay))
                            visited.add((ax, ay))
                        elif board[ax][ay] == '.':
                            liberties.add((ax, ay))
            
            # If no liberties, the group is captured
            if not liberties:
                captured.update(visited)
    
    # Remove temporary stone
    board[x][y] = '.'
    
    return len(captured)

def find_best_move():
    board = create_board()
    best_move = None
    max_captures = 0
    
    # Try each empty position
    for x in range(11):
        for y in range(11):
            if board[x][y] == '.':
                captures = count_captures(board, x, y)
                if captures > max_captures:
                    max_captures = captures
                    best_move = (x, y)
    
    if best_move:
        # Convert to board coordinates (A-K, 1-11)
        col = chr(best_move[1] + ord('A'))
        row = 11 - best_move[0]
        return f"{col}{row}", max_captures
    
    return None, 0

move, captures = find_best_move()
print(f"Best move: {move}")
print(f"Captures: {captures}")