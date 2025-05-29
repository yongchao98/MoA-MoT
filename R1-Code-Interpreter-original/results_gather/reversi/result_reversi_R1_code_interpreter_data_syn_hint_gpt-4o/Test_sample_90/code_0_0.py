def print_board(board):
    for row in board:
        print(''.join(row))

def is_on_board(x, y, n):
    return 0 <= x < n and 0 <= y < n

def flip_pieces(board, start_x, start_y, dx, dy, player):
    x, y = start_x, start_y
    pieces_to_flip = []
    while True:
        x += dx
        y += dy
        if not is_on_board(x, y, len(board)) or board[x][y] == '*':
            return
        if board[x][y] == player:
            for px, py in pieces_to_flip:
                board[px][py] = player
            return
        pieces_to_flip.append((x, y))

def make_move(board, x, y, player):
    board[x][y] = player
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (-1, 1), (1, -1), (1, 1)]
    for dx, dy in directions:
        flip_pieces(board, x, y, dx, dy, player)

# Initial board setup
board = [
    ['*', '*', '*', '*', '*', '*'],
    ['*', '*', '*', '*', '*', '*'],
    ['*', '*', '1', '0', '*', '*'],
    ['*', '*', '0', '1', '0', '*'],
    ['*', '*', '*', '*', '0', '*'],
    ['*', '*', '*', '*', '0', '*']
]

# Moves
moves = [(3, 4, '0'), (4, 4, '1'), (5, 4, '0')]

# Apply moves
for x, y, player in moves:
    make_move(board, x - 1, y - 1, player)

# Print final board
print_board(board)