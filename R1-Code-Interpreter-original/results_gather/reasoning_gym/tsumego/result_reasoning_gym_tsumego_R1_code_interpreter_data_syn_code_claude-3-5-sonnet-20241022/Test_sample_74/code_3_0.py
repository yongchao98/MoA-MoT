def print_board_with_move(move_x, move_y):
    # Initialize board
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
    
    # Place the move
    board[move_x][move_y] = 'M'
    
    # Print board with coordinates
    print("   A B C D E F G H I J K L")
    for i in range(12):
        row_num = str(12-i).rjust(2)
        row = ' '.join(board[i])
        print(f"{row_num} {row}")
    
    # Print surrounding stones for the move
    move_col = chr(move_y + ord('A'))
    move_row = 12 - move_x
    print(f"\nAnalyzing move at {move_col}{move_row}")

# Test H4 position
print("Testing H4 position:")
print_board_with_move(8, 7)