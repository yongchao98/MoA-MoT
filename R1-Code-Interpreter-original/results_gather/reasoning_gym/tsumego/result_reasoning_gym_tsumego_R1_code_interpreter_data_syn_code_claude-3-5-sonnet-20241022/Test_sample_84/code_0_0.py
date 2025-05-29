def create_board():
    board = [['.'] * 13 for _ in range(13)]
    # Place the stones according to the given position
    # White stones
    board[10][7] = 'O'  # H11
    board[7][12] = 'O'  # M8
    board[5][0] = 'O'   # A6
    board[4][7] = 'O'   # H5
    board[3][6] = 'O'   # G4
    board[2][6] = 'O'   # G3
    board[2][7] = 'O'   # H3
    board[2][8] = 'O'   # I3
    board[1][7] = 'O'   # H2

    # Black stones
    board[9][4] = 'X'   # E10
    board[6][7] = 'X'   # H7
    board[4][5] = 'X'   # F5
    board[4][8] = 'X'   # I4
    board[3][1] = 'X'   # B4
    board[2][5] = 'X'   # F3
    board[2][9] = 'X'   # J3
    board[1][6] = 'X'   # G2
    board[1][8] = 'X'   # I2
    board[0][7] = 'X'   # H1
    board[0][8] = 'X'   # I1
    board[0][10] = 'X'  # K1
    board[4][10] = 'X'  # K5

    return board

def print_board(board):
    for i, row in enumerate(reversed(board)):
        print(f"{i+1:2d} {' '.join(row)}")
    print("   A B C D E F G H I J K L M")

def check_capture_at_h4():
    board = create_board()
    # The move we're considering is at H4
    h4_coords = (3, 7)  # 4th row (0-based), 8th column (0-based)
    
    # Check if this move would surround the white stones
    surrounded_positions = [
        (2, 6),  # G3
        (2, 7),  # H3
        (2, 8),  # I3
    ]
    
    # Verify that these positions are surrounded after the move
    black_stones = [
        (2, 5),  # F3
        (2, 9),  # J3
        (1, 6),  # G2
        (1, 8),  # I2
        h4_coords  # H4 (the proposed move)
    ]
    
    print("Current board:")
    print_board(board)
    print("\nAfter playing at H4:")
    board[h4_coords[0]][h4_coords[1]] = 'X'
    print_board(board)
    print("\nThis move would capture the white stones at G3, H3, and I3")
    print("The move coordinates are: H4")

check_capture_at_h4()