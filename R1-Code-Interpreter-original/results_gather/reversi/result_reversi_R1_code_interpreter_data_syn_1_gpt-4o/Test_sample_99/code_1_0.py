def print_board(board):
    for row in board:
        print(''.join(row))

def is_on_board(x, y, n):
    return 0 <= x < n and 0 <= y < n

def flip_pieces(board, x, y, player, opponent, n):
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (-1, 1), (1, -1), (1, 1)]
    to_flip = []
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        pieces_to_flip = []
        while is_on_board(nx, ny, n) and board[nx][ny] == opponent:
            pieces_to_flip.append((nx, ny))
            nx += dx
            ny += dy
        if is_on_board(nx, ny, n) and board[nx][ny] == player:
            to_flip.extend(pieces_to_flip)
    for fx, fy in to_flip:
        board[fx][fy] = player

def apply_move(board, move, player, opponent, n):
    x, y = move
    board[x][y] = player
    flip_pieces(board, x, y, player, opponent, n)

# Initial board setup
n = 6
board = [
    ['*', '*', '*', '*', '*', '*'],
    ['*', '1', '0', '*', '*', '*'],
    ['*', '0', '0', '0', '*', '*'],
    ['*', '*', '0', '1', '*', '*'],
    ['*', '*', '*', '*', '*', '*'],
    ['*', '*', '*', '*', '*', '*']
]

# Moves
moves = [(2, 3), (2, 2), (3, 2)]
players = ['0', '1', '0']  # Black, White, Black

# Apply moves
for i, move in enumerate(moves):
    player = players[i]
    opponent = '1' if player == '0' else '0'
    apply_move(board, move, player, opponent, n)

# Print final board state
print_board(board)