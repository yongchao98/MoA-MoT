def check_liberties(board, x, y, color, visited=None):
    if visited is None:
        visited = set()
    
    if (x, y) in visited or x < 0 or y < 0 or x >= 13 or y >= 13:
        return set()
    if board[y][x] != color:
        return set()
    
    visited.add((x, y))
    liberties = set()
    
    # Check all adjacent positions
    for dx, dy in [(-1,0), (1,0), (0,-1), (0,1)]:
        new_x, new_y = x + dx, y + dy
        if new_x < 0 or new_y < 0 or new_x >= 13 or new_y >= 13:
            continue
        
        if board[new_y][new_x] == '.':
            liberties.add((new_x, new_y))
        elif board[new_y][new_x] == color:
            liberties.update(check_liberties(board, new_x, new_y, color, visited))
    
    return liberties

# Create the board
board = [['.'] * 13 for _ in range(13)]

# Set up the relevant stones
# White stones
board[10][0] = 'O'  # A11
board[7][2] = 'O'   # C8
board[7][5] = 'O'   # F8
board[7][6] = 'O'   # G8
board[7][7] = 'O'   # H8
board[6][6] = 'O'   # G7
board[4][3] = 'O'   # D5
board[4][4] = 'O'   # E5
board[4][8] = 'O'   # I5

# Black stones
board[11][6] = 'X'  # G12
board[10][7] = 'X'  # H11
board[8][5] = 'X'   # F9
board[8][7] = 'X'   # H9
board[7][4] = 'X'   # E8
board[7][8] = 'X'   # I8
board[6][2] = 'X'   # C7
board[6][5] = 'X'   # F7
board[6][7] = 'X'   # H7
board[5][6] = 'X'   # G6
board[2][3] = 'X'   # D3
board[0][10] = 'X'  # K1

# Test placing black stone at G6
board[5][6] = 'X'   # Place black stone at G6

# Check liberties of the white group
white_group_liberties = check_liberties(board, 6, 7, 'O')  # Check from G7
print(f"Liberties for white group after Black plays at G6: {white_group_liberties}")
if not white_group_liberties:
    print("The white group would be captured!")