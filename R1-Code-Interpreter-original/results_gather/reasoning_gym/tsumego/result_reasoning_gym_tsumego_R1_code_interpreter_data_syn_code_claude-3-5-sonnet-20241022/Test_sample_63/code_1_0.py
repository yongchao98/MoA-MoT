def create_board():
    board = [['.'] * 13 for _ in range(13)]
    # Place White stones (O) - adjusting coordinates to 0-based indexing from bottom-left
    white_positions = [
        (0, 12), # A1
        (3, 12), # D1
        (2, 11), # C2
        (3, 11), # D2
        (5, 11), # F2
        (3, 10), # D3
        (4, 10), # E3
        (3, 9),  # D4
        (11, 8), # L5
        (9, 0)   # J13
    ]
    
    # Place Black stones (X)
    black_positions = [
        (2, 12), # C1
        (4, 12), # E1
        (6, 12), # G1
        (1, 11), # B2
        (2, 10), # C3
        (5, 1),  # F12
        (2, 2),  # C11
        (7, 2),  # H11
        (9, 5),  # J8
        (9, 6)   # J7
    ]
    
    for x, y in white_positions:
        board[y][x] = 'O'
    for x, y in black_positions:
        board[y][x] = 'X'
    return board

def is_valid(x, y):
    return 0 <= x < 13 and 0 <= y < 13

def count_captures(board, x, y):
    if board[y][x] != '.':
        return 0
    
    directions = [(0,1), (1,0), (0,-1), (-1,0)]
    captured = set()
    
    # Temporarily place the stone
    board[y][x] = 'X'
    
    for dx, dy in directions:
        stones = []
        nx, ny = x + dx, y + dy
        
        while is_valid(nx, ny) and board[ny][nx] == 'O':
            stones.append((nx, ny))
            nx += dx
            ny += dy
            
        if stones and is_valid(nx, ny) and board[ny][nx] == 'X':
            captured.update(stones)
    
    # Remove the temporary stone
    board[y][x] = '.'
    return len(captured)

def find_best_move():
    board = create_board()
    best_move = None
    max_captures = 0
    
    for y in range(13):
        for x in range(13):
            if board[y][x] == '.':
                captures = count_captures(board, x, y)
                if captures > max_captures:
                    max_captures = captures
                    best_move = (x, y)
    
    if best_move:
        # Convert to board coordinates (A-M, 1-13)
        col = chr(ord('A') + best_move[0])
        row = 13 - best_move[1]
        print(f"Best move: {col}{row} with {max_captures} captures")

find_best_move()