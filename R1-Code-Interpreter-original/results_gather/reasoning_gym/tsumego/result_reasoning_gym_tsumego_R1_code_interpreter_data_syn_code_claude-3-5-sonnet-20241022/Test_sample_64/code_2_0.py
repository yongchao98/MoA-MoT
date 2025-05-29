def analyze_vital_points():
    # Create board representation
    board = [['.'] * 11 for _ in range(11)]
    
    # Set up position
    black = [(0,0), (1,0), (6,8), (7,8), (5,8), (5,6), (7,6), (5,3), (5,1), (8,4), (10,6)]
    white = [(1,7), (3,6), (6,7), (6,8), (6,6), (8,7), (8,6), (0,3)]
    
    for x, y in black:
        board[10-y][x] = 'X'
    for x, y in white:
        board[10-y][x] = 'O'
    
    # Check potential cutting points and vital points
    vital_points = []
    for x in range(11):
        for y in range(11):
            if board[y][x] == '.':
                # Count adjacent stones
                black_adj = 0
                white_adj = 0
                for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
                    nx, ny = x + dx, y + dy
                    if 0 <= nx < 11 and 0 <= ny < 11:
                        if board[ny][nx] == 'X':
                            black_adj += 1
                        elif board[ny][nx] == 'O':
                            white_adj += 1
                
                if black_adj >= 1 and white_adj >= 1:  # Potential vital point
                    vital_points.append((x, y))
    
    print("Vital points analysis:")
    for x, y in vital_points:
        print(f"({chr(65+x)},{11-y}): Adjacent to {board[y-1][x] if y>0 else '.'} (N), "
              f"{board[y+1][x] if y<10 else '.'} (S), "
              f"{board[y][x-1] if x>0 else '.'} (W), "
              f"{board[y][x+1] if x<10 else '.'} (E)")

analyze_vital_points()