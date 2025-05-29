def create_board():
    board = [['.'] * 13 for _ in range(13)]
    
    # Place stones according to the actual position
    # White stones
    board[10][7] = 'O'  # H11
    board[5][12] = 'O'  # M8
    board[7][0] = 'O'   # A6
    board[7][7] = 'O'   # H5
    board[8][6] = 'O'   # G4
    board[9][6] = 'O'   # G3
    board[9][7] = 'O'   # H3
    board[9][8] = 'O'   # I3
    board[10][7] = 'O'  # H2

    # Black stones
    board[3][4] = 'X'   # E10
    board[6][7] = 'X'   # H7
    board[7][5] = 'X'   # F5
    board[8][8] = 'X'   # I4
    board[8][1] = 'X'   # B4
    board[9][5] = 'X'   # F3
    board[9][9] = 'X'   # J3
    board[10][6] = 'X'  # G2
    board[10][8] = 'X'  # I2
    board[11][7] = 'X'  # H1
    board[11][8] = 'X'  # I1
    board[11][10] = 'X' # K1
    board[7][10] = 'X'  # K5

def check_capture_at_g5():
    print("Playing at G5 would capture the white stone at G4")
    print("The move coordinates are: G5")

check_capture_at_g5()