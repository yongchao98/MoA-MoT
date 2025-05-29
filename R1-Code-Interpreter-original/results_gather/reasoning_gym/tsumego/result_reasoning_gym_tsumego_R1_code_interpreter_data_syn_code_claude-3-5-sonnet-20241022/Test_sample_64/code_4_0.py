def analyze_position():
    # Create board representation
    board = [['.'] * 11 for _ in range(11)]
    
    # Set up position
    black = [(0,0), (1,0), (6,8), (7,8), (5,8), (5,6), (7,6), (5,3), (5,1), (8,4), (10,6)]
    white = [(1,7), (3,6), (6,7), (6,8), (6,6), (8,7), (8,6), (0,3)]
    
    for x, y in black:
        board[10-y][x] = 'X'
    for x, y in white:
        board[10-y][x] = 'O'
    
    # Key points to analyze
    key_points = [(5,7), (6,5), (7,7)]  # F7, G5, H7
    
    for x, y in key_points:
        # Make a copy of the board
        test_board = [row[:] for row in board]
        test_board[10-y][x] = 'X'
        
        # Check if this move threatens the white group
        white_liberties = 0
        white_connections = 0
        
        # Check liberties and connections for critical white stones
        for wx, wy in [(6,7), (6,6)]:  # G7, G6
            for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
                nx, ny = wx + dx, wy + dy
                if 0 <= nx < 11 and 0 <= ny < 11:
                    if test_board[10-ny][nx] == '.':
                        white_liberties += 1
                    elif test_board[10-ny][nx] == 'O':
                        white_connections += 1
        
        print(f"Move at ({chr(65+x)},{y})")
        print(f"White liberties: {white_liberties}")
        print(f"White connections: {white_connections}")
        print("---")

analyze_position()