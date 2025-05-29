def print_board(board):
    print("   A B C D E F G H I J K L M")
    for i in range(13, 0, -1):
        print(f"{i:2d}", end=" ")
        for j in range(13):
            print(board[i-1][j], end=" ")
        print()

def check_liberties(x, y, board):
    directions = [(0,1), (0,-1), (1,0), (-1,0)]
    liberties = 0
    for dx, dy in directions:
        new_x, new_y = x + dx, y + dy
        if 0 <= new_x < 13 and 0 <= new_y < 13:
            if board[new_y][new_x] == '.':
                liberties += 1
    return liberties

# Create board representation
board = [
    ['X', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', 'X', '.', '.', '.', '.', 'O', '.', '.', '.', '.'],
    ['.', '.', 'X', 'O', 'X', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', 'X', 'O', 'O', 'O', 'X', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', 'X', '.', 'X', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', 'O', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', 'X', '.', '.', '.', '.', '.', '.', '.', 'X', '.', '.', '.'],
    ['.', '.', '.', '.', 'X', 'X', '.', '.', '.', '.', '.', '.', '.'],
    ['O', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'O', '.', '.', 'O'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', 'X']
]

# Check the critical white group
white_group = [(2,3), (3,3), (4,3), (3,2)]  # D3, D4, E4, C4 coordinates
empty_adjacent = set()

for x, y in white_group:
    directions = [(0,1), (0,-1), (1,0), (-1,0)]
    for dx, dy in directions:
        new_x, new_y = x + dx, y + dy
        if 0 <= new_x < 13 and 0 <= new_y < 13:
            if board[new_y][new_x] == '.':
                empty_adjacent.add((new_x, new_y))

print("Critical empty points adjacent to white group:")
for x, y in empty_adjacent:
    col = chr(65 + x)
    row = y + 1
    print(f"{col}{row}")