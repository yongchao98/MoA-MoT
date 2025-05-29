def print_board(board):
    print("   A B C D E F G H I J K L")
    for i in range(12, 0, -1):
        row = f"{i:2d} "
        for j in range(12):
            row += board[i-1][j] + " "
        print(row)

def create_board():
    board = [['.'] * 12 for _ in range(12)]
    # Place black stones (X)
    black_stones = [(9,8), (9,9), (9,10), (9,11),
                    (7,8), (6,9), (6,10),
                    (5,1), (5,2),
                    (0,0), (0,5)]
    # Place white stones (O)
    white_stones = [(11,8),
                    (8,9), (8,10), (8,11),
                    (7,10),
                    (5,4)]
    
    for x, y in black_stones:
        board[x][y] = 'X'
    for x, y in white_stones:
        board[x][y] = 'O'
    
    return board

# Create and print the board
board = create_board()
print_board(board)

# Key points to consider for capturing move
key_points = [(8,11), (7,11), (7,9)]  # L9, L8, J8 coordinates
print("\nKey points to consider:")
for x, y in key_points:
    col = chr(y + ord('A'))
    row = x + 1
    print(f"{col}{row}")