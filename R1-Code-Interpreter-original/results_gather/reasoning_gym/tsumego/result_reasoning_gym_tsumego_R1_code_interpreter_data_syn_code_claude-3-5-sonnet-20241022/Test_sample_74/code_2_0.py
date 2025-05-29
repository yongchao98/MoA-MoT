def check_specific_moves():
    # Initialize board with key stones around H4 and H7 areas
    board = [['.'] * 12 for _ in range(12)]
    
    # Place stones exactly as in the problem
    # Black stones (X)
    black = [
        (11,0),  # A11
        (10,4), (10,6), (10,7),  # E10, G10, H10
        (9,5), (9,8), (9,11),  # F9, I9, L9
        (8,1), (8,5),  # B8, F8
        (7,1), (7,2), (7,6),  # B7, C7, G7
        (5,5),  # F5
        (1,10)  # K1
    ]
    
    # White stones (O)
    white = [
        (11,8),  # I12
        (9,4), (9,6), (9,7),  # E9, G9, H9
        (8,6),  # G8
        (5,0),  # A5
        (1,7)   # H1
    ]
    
    for x, y in black:
        board[x][y] = 'X'
    for x, y in white:
        board[x][y] = 'O'
    
    def print_local_area(x, y):
        print(f"Area around {chr(y+ord('A'))}{12-x}:")
        for i in range(max(0,x-2), min(12,x+3)):
            row = ''
            for j in range(max(0,y-2), min(12,y+3)):
                row += board[i][j] + ' '
            print(row)
    
    # Check H7 (7,7)
    print_local_area(7,7)
    
    # Check H4 (8,7)
    print_local_area(8,7)

print("Analyzing specific positions:")
check_specific_moves()