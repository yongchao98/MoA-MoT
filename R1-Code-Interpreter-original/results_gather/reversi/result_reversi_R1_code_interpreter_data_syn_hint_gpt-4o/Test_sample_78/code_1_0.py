# Initial board setup
board = [
    ['*', '*', '*', '*', '*', '*'],
    ['*', '*', '*', '*', '*', '*'],
    ['*', '*', '1', '1', '1', '*'],
    ['*', '*', '0', '0', '0', '*'],
    ['*', '*', '*', '*', '*', '*'],
    ['*', '*', '*', '*', '*', '*']
]

# Directions for checking: right, left, down, up, diagonal down-right, diagonal up-left, diagonal down-left, diagonal up-right
directions = [(0, 1), (0, -1), (1, 0), (-1, 0), (1, 1), (-1, -1), (1, -1), (-1, 1)]

def is_valid(x, y):
    return 0 <= x < 6 and 0 <= y < 6

def flip_pieces(board, x, y, player):
    opponent = '1' if player == '0' else '0'
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        pieces_to_flip = []
        while is_valid(nx, ny) and board[nx][ny] == opponent:
            pieces_to_flip.append((nx, ny))
            nx += dx
            ny += dy
        if is_valid(nx, ny) and board[nx][ny] == player:
            for px, py in pieces_to_flip:
                board[px][py] = player

# Round 1: Black's move at (4,5)
board[3][4] = '0'
flip_pieces(board, 3, 4, '0')

# Round 2: White's move at (3,5)
board[2][4] = '1'
flip_pieces(board, 2, 4, '1')

# Flatten the board to a single list
flattened_board = [cell for row in board for cell in row]
print(flattened_board)