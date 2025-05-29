def check_capture(move):
    # Board representation (0-based indexing)
    board = [
        ['.','.','.','.','.','.','.','.','.','.','.'],
        ['.','.','.','.','.','.','.','.','.','.','.'],
        ['.','.','.','.','.','.','.','.','.','.','O'],
        ['.','.','.','.','.','.','.','.','.','X','.'],
        ['.','.','X','.','.','.','.','.','.','.','X'],
        ['.','.','.','.','.','.','.','.','X','O','O'],
        ['.','.','.','.','.','.','.','.','X','O','.'],
        ['.','.','.','.','.','.','.','.','.','.','X'],
        ['.','.','O','.','.','.','.','.','.','.','O'],
        ['.','O','.','.','X','.','.','.','.','.','.'],
        ['.','.','.','.','.','X','.','.','.','.','.']
    ]
    
    x, y = move
    captures = 0
    
    # Check if move at I7 captures the J7-K7 group
    if board[x][y] == '.' and x == 7 and y == 8:  # I7 in zero-based coordinates
        # Check if surrounded
        if (board[x][y+1] == 'O' and board[x][y+2] == 'O' and  # J7-K7
            board[x-1][y+1] == 'X' and  # J8
            board[x+1][y+1] == 'X'):    # J6
            captures = 2
    
    return captures

# Test the I7 move (coordinates in 0-based indexing)
move = (7, 8)  # I7 in 0-based coordinates
print(f"Stones captured by playing at I7: {check_capture(move)}")