def print_board(board):
    print("   A B C D E F G H I J K L M")
    for i in range(13, 0, -1):
        print(f"{i:2d}", end=" ")
        for j in range(13):
            print(board[i-1][j], end=" ")
        print()

def get_liberties(board, x, y):
    liberties = 0
    directions = [(0,1), (0,-1), (1,0), (-1,0)]
    for dx, dy in directions:
        new_x, new_y = x + dx, y + dy
        if 0 <= new_x < 13 and 0 <= new_y < 13:
            if board[new_y][new_x] == '.':
                liberties += 1
    return liberties

# Initialize the board
board = [
    ['.', '.', '.', '.', '.', 'X', 'O', 'X', '.', '.', '.', '.', '.'],  # 1
    ['.', '.', '.', '.', 'X', 'O', 'O', 'O', 'X', '.', '.', '.', '.'],  # 2
    ['.', '.', '.', '.', '.', 'X', '.', 'X', '.', '.', '.', '.', '.'],  # 3
    ['.', '.', '.', '.', '.', '.', 'X', '.', '.', '.', '.', '.', '.'],  # 4
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],  # 5
    ['O', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],  # 6
    ['.', '.', '.', '.', '.', '.', 'X', 'X', '.', '.', '.', '.', '.'],  # 7
    ['.', '.', '.', 'X', '.', '.', '.', '.', '.', '.', '.', '.', '.'],  # 8
    ['.', '.', '.', '.', '.', '.', '.', 'O', '.', '.', '.', '.', '.'],  # 9
    ['.', '.', '.', '.', '.', 'O', '.', '.', 'O', 'X', '.', '.', '.'],  # 10
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],  # 11
    ['.', '.', '.', '.', 'X', '.', '.', 'O', '.', '.', 'O', '.', '.'],  # 12
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.']   # 13
]

# Check liberties of the white stones at F2-H2
white_group = [(5,1), (6,1), (7,1)]  # F2, G2, H2 in zero-based coordinates
total_liberties = 0
for x, y in white_group:
    liberties = get_liberties(board, x, y)
    total_liberties += liberties
    print(f"White stone at ({chr(65+x)}{y+1}) has {liberties} liberties")

# Check if G1 is the key point
g1_x, g1_y = 6, 0  # G1 in zero-based coordinates
if board[g1_y][g1_x] == 'O':
    g1_liberties = get_liberties(board, g1_x, g1_y)
    print(f"\nWhite stone at G1 has {g1_liberties} liberties")