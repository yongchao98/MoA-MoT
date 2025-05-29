def print_board_with_move(move_x, move_y):
    # Initialize board with correct coordinates (12 is top row)
    board = [['.'] * 12 for _ in range(12)]
    
    # Place black stones (X)
    black = [
        (11,0),  # A11
        (10,4), (10,6), (10,7),  # E10, G10, H10
        (9,5), (9,8), (9,11),  # F9, I9, L9
        (8,1), (8,5),  # B8, F8
        (7,1), (7,2), (7,6),  # B7, C7, G7
        (5,5),  # F5
        (1,10)  # K1
    ]
    
    # Place white stones (O)
    white = [
        (11,8),  # I12
        (9,4), (9,6), (9,7),  # E9, G9, H9
        (8,6),  # G8
        (5,0),  # A5
        (1,7)   # H1
    ]
    
    for x, y in black:
        board[12-x][y] = 'X'
    for x, y in white:
        board[12-x][y] = 'O'
    
    # Mark the move
    board[12-move_x][move_y] = 'M'
    
    # Print board
    print("   A B C D E F G H I J K L")
    for i in range(12):
        row_num = str(12-i).rjust(2)
        row = ' '.join(board[i])
        print(f"{row_num} {row}")
    
    # Print analysis of surrounding area
    move_col = chr(move_y + ord('A'))
    print(f"\nAnalyzing move at {move_col}{move_x}")

# Test F7 position
print("Testing F7 position:")
print_board_with_move(7, 5)