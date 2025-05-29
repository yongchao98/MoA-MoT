def create_board():
    board = {}
    for i in range(1, 12):
        for c in 'ABCDEFGHIJK':
            board[f"{c}{i}"] = '.'
    
    # Set white stones (O)
    white_stones = ['A11', 'H9', 'G6', 'H5', 'I5', 'D4', 'H4', 'I4', 'I3']
    for pos in white_stones:
        board[pos] = 'O'
    
    # Set black stones (X)
    black_stones = ['B7', 'I7', 'I6', 'C5', 'J5', 'G4', 'J4', 'H3', 'J3', 'I2', 'C1', 'G1']
    for pos in black_stones:
        board[pos] = 'X'
    
    return board

def get_liberties(board, pos):
    col, row = pos[0], int(pos[1:])
    liberties = []
    # Check all adjacent positions
    for dc, dr in [(0,1), (0,-1), (1,0), (-1,0)]:
        new_col = chr(ord(col) + dc)
        new_row = row + dr
        if 'A' <= new_col <= 'K' and 1 <= new_row <= 11:
            new_pos = f"{new_col}{new_row}"
            if new_pos in board and board[new_pos] == '.':
                liberties.append(new_pos)
    return liberties

board = create_board()
# Check each white stone and its liberties
for pos in ['G6', 'H5', 'I5', 'H4', 'I4', 'I3']:
    libs = get_liberties(board, pos)
    print(f"White stone at {pos} has liberties at: {libs}")