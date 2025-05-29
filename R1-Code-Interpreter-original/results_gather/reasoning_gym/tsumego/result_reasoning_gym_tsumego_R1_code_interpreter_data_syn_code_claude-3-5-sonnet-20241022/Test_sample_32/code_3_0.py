def create_board():
    board = [['.'] * 12 for _ in range(12)]
    
    # Place stones using exact coordinates from the diagram
    # Black stones (X)
    black = [
        (11,1),  # B12
        (9,9), (9,10),  # J10, K10
        (8,8), (8,11),  # I9, L9
        (7,8),  # I8
        (5,3), (5,7),  # D6, H6
        (3,10),  # K4
        (2,2)   # C2
    ]
    
    # White stones (O)
    white = [
        (10,0), (10,10),  # A11, K11
        (8,9), (8,10),  # J9, K9
        (7,11),  # L8
        (7,3),  # D8
        (2,0),  # A3
        (1,3)   # D2
    ]
    
    for x, y in black:
        board[x][y] = 'X'
    for x, y in white:
        board[x][y] = 'O'
    
    # Analyze K8 position (7,10)
    captures = 0
    x, y = 7, 10  # K8 position
    
    # Check if placing a stone here would capture any groups
    adjacent = [(x+1,y), (x-1,y), (x,y+1), (x,y-1)]
    white_groups = []
    
    for ax, ay in adjacent:
        if (0 <= ax < 12 and 0 <= ay < 12 and 
            board[ax][ay] == 'O'):
            white_groups.append((ax,ay))
    
    print(f"Position K8 is adjacent to {len(white_groups)} white stones")
    print(f"White stones at: {[(12-x, chr(ord('A')+y)) for x,y in white_groups]}")

create_board()