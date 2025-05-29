def print_board(board):
    print("   A B C D E F G H I J K L M")
    for i in range(13, 0, -1):
        print(f"{i:2d}", end=" ")
        for j in range(13):
            print(board[i-1][j], end=" ")
        print()

def check_move(board, x, y, move):
    # Copy board and make move
    new_board = [row[:] for row in board]
    if move == 'B2':
        new_board[1][1] = 'X'
    elif move == 'C2':
        new_board[1][2] = 'X'
    
    # Check surrounding liberties
    liberties = 0
    adjacent = [(0,1), (0,-1), (1,0), (-1,0)]
    for dx, dy in adjacent:
        nx, ny = x + dx, y + dy
        if 0 <= nx < 13 and 0 <= ny < 13 and new_board[ny][nx] == '.':
            liberties += 1
    return liberties

# Initialize the board
board = [['.'] * 13 for _ in range(13)]

# Set up the relevant stones
stones = {
    'X': [(0,1), (1,2), (2,3), (3,2), (3,0)],  # Black stones
    'O': [(2,0), (2,1), (2,2), (4,2)]   # White stones
}

for color in stones:
    for x, y in stones[color]:
        board[y][x] = color

# Check potential moves
moves = ['B2']
for move in moves:
    print(f"\nAnalyzing move {move}:")
    liberties = check_move(board, 1, 1, move)
    print(f"Resulting liberties after {move}: {liberties}")