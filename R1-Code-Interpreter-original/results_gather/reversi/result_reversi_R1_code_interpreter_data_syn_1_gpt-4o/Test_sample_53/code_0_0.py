# Initial board configuration
board = [
    ['*', '*', '*', '*'],
    ['*', '1', '0', '*'],
    ['*', '0', '0', '*'],
    ['0', '0', '0', '*']
]

# Function to flip pieces in a given direction
def flip_pieces(board, start, direction, player):
    x, y = start
    dx, dy = direction
    to_flip = []
    x += dx
    y += dy
    while 0 <= x < 4 and 0 <= y < 4:
        if board[x][y] == '*':
            return
        if board[x][y] == player:
            for fx, fy in to_flip:
                board[fx][fy] = player
            return
        to_flip.append((x, y))
        x += dx
        y += dy

# Function to place a piece and flip opponent's pieces
def place_piece(board, position, player):
    x, y = position
    board[x][y] = player
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (1, 1), (-1, 1), (1, -1)]
    for direction in directions:
        flip_pieces(board, position, direction, player)

# Round 1: Your move (Black) at (4,3)
place_piece(board, (3, 2), '0')

# Round 2: Opponent's move (White) at (4,2)
place_piece(board, (3, 1), '1')

# Round 3: Your move (Black) at (4,1)
place_piece(board, (3, 0), '0')

# Flatten the board to a single list for the final output
flattened_board = [cell for row in board for cell in row]
print(flattened_board)