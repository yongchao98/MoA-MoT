def analyze_move(move_x, move_y):
    # Create board representation
    board = [['.'] * 11 for _ in range(11)]
    
    # Set up position
    black = [(0,0), (1,0), (6,8), (7,8), (5,8), (5,6), (7,6), (5,3), (5,1), (8,4), (10,6)]
    white = [(1,7), (3,6), (6,7), (6,8), (6,6), (8,7), (8,6), (0,3)]
    
    for x, y in black:
        board[10-y][x] = 'X'
    for x, y in white:
        board[10-y][x] = 'O'
    
    # Place the move
    board[10-move_y][move_x] = 'X'
    
    # Check impact on white group
    white_group_liberties = set()
    white_core = [(6,7), (6,8)]  # Focus on these critical white stones
    
    for wx, wy in white_core:
        for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
            nx, ny = wx + dx, wy + dy
            if 0 <= nx < 11 and 0 <= ny < 11 and board[ny][nx] == '.':
                white_group_liberties.add((nx, ny))
    
    print(f"Move at ({chr(65+move_x)},{move_y})")
    print(f"Critical white stones liberties: {len(white_group_liberties)}")
    print(f"Liberty points: {[(chr(65+x), 11-y) for x,y in white_group_liberties]}")

analyze_move(5, 8)  # Analyzing F8