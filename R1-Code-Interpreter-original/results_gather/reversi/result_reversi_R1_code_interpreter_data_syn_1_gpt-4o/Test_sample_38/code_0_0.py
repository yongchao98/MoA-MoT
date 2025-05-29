# Initial board setup
board = [
    ['*', '*', '*', '*'],
    ['*', '1', '1', '1'],
    ['*', '0', '0', '*'],
    ['*', '*', '0', '*']
]

# Function to flip pieces in a given direction
def flip_pieces(board, start, direction, player):
    x, y = start
    dx, dy = direction
    to_flip = []
    opponent = '0' if player == '1' else '1'
    
    while 0 <= x < 4 and 0 <= y < 4:
        x += dx
        y += dy
        if 0 <= x < 4 and 0 <= y < 4:
            if board[x][y] == opponent:
                to_flip.append((x, y))
            elif board[x][y] == player:
                for fx, fy in to_flip:
                    board[fx][fy] = player
                break
            else:
                break

# Function to make a move
def make_move(board, position, player):
    x, y = position
    board[x][y] = player
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (1, 1), (-1, 1), (1, -1)]
    for direction in directions:
        flip_pieces(board, (x, y), direction, player)

# Round 1: Black's move
make_move(board, (3, 2), '0')

# Round 2: White's move
make_move(board, (1, 3), '1')

# Flatten the board to a single list
flattened_board = [cell for row in board for cell in row]
print(flattened_board)