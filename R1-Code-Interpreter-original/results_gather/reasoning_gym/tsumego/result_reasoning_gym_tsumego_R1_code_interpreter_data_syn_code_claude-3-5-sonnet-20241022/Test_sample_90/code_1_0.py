def print_board(board):
    print("   A B C D E F G H I J K L M")
    for i in range(13):
        row = 13 - i
        print(f"{row:2d} {' '.join(board[i])}")

def get_liberties(board, x, y):
    liberties = set()
    directions = [(0,1), (0,-1), (1,0), (-1,0)]
    for dx, dy in directions:
        new_x, new_y = x + dx, y + dy
        if 0 <= new_x < 13 and 0 <= new_y < 13:
            if board[new_y][new_x] == '.':
                liberties.add((new_x, new_y))
    return liberties

# Create the board
board = [['.'] * 13 for _ in range(13)]

# Place the stones as per the diagram
board[12-12][11] = 'X'  # M1
board[12-11][10] = 'X'  # L2
board[12-3][10] = 'O'   # K3
board[12-7][11] = 'O'   # L7
board[12-7][12] = 'O'   # M7
board[12-10][7] = 'O'   # H10
board[12-10][9] = 'O'   # J10
board[12-6][5] = 'X'    # F6
board[12-6][6] = 'X'    # G6
board[12-5][3] = 'X'    # D5
board[12-5][4] = 'X'    # E5
board[12-5][5] = 'O'    # F5
board[12-5][6] = 'O'    # G5
board[12-5][7] = 'X'    # H5
board[12-4][4] = 'X'    # E4
board[12-4][5] = 'O'    # F4
board[12-3][5] = 'X'    # F3
board[12-2][3] = 'O'    # D2

# Check liberties of the white group
white_group = [(5,8), (5,7), (6,7)]  # F5, F4, G5
all_liberties = set()
for x, y in white_group:
    libs = get_liberties(board, x, y)
    all_liberties.update(libs)
    print(f"Liberties at ({x},{y}): {libs}")

print("\nTotal group liberties:", all_liberties)
print("\nBoard position:")
print_board(board)