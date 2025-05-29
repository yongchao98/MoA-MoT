def print_board(move=None):
    board = [
        ['.'] * 13 for _ in range(13)
    ]
    
    # Place existing stones accurately
    black_stones = [
        (0,7), (1,6), (0,5), (3,6), # left side
        (1,4), (2,3), (5,3), # middle area
        (10,2), (11,2), # right side
        (5,0)  # bottom
    ]
    
    white_stones = [
        (2,7), (4,6), # center group
        (1,5), (2,5), (3,5), # center group
        (7,7), (9,7), # top stones
        (7,2), (10,2) # bottom stones
    ]
    
    for x,y in black_stones:
        board[12-y][x] = 'X'
    for x,y in white_stones:
        board[12-y][x] = 'O'
        
    if move:
        x, y = move
        board[12-y][x] = '*'
    
    print('   A B C D E F G H I J K L M')
    for i in range(13):
        print(f"{13-i:2d} {' '.join(board[i])}")

def analyze_position():
    # The central white group (B6-D6, C5) is the main target
    # Looking at E7 as the key move
    print("\nAnalyzing key move at E7:")
    print_board((4,6))  # E7 coordinates
    
    print("\nThis move would capture the following white stones:")
    print("- B6 (O)")
    print("- C6 (O)")
    print("- D6 (O)")
    print("- C5 (O)")
    print("\nTotal captures: 4 stones")

analyze_position()