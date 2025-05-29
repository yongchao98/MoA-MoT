def print_board(board):
    print("   A B C D E F G H I J K L M")
    for i in range(13):
        row = 13 - i
        print(f"{row:2d} {' '.join(board[i])}")

def get_liberties(board, x, y):
    liberties = 0
    directions = [(0,1), (0,-1), (1,0), (-1,0)]
    for dx, dy in directions:
        new_x, new_y = x + dx, y + dy
        if 0 <= new_x < 13 and 0 <= new_y < 13:
            if board[new_y][new_x] == '.':
                liberties += 1
    return liberties

# Create the board
board = [['.'] * 13 for _ in range(13)]

# Place the stones as per the diagram
# Black stones (X)
stones_black = [(11,1), (10,1), (5,3), (4,3), (5,4), (4,4), (6,5), (6,6), (5,7), 
                (10,11), (12,11), (11,12), (0,12)]

# White stones (O)
stones_white = [(7,9), (9,9), (6,4), (5,4), (6,5), (6,6), (11,6), (12,6), (2,1)]

for x, y in stones_black:
    board[y][x] = 'X'
for x, y in stones_white:
    board[y][x] = 'O'

# Check liberties of white stones at F4-F5-G5
critical_stones = [(5,4), (5,5), (6,5)]
for x, y in critical_stones:
    liberties = get_liberties(board, x, y)
    print(f"Liberties at ({x},{y}): {liberties}")

print_board(board)