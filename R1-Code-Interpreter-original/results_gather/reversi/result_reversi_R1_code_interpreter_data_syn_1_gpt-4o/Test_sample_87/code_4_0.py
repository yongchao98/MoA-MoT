# Initial board configuration
board = [
    ['*', '*', '*', '*', '*', '*'],
    ['*', '*', '0', '*', '*', '*'],
    ['*', '*', '0', '0', '*', '*'],
    ['*', '*', '0', '1', '*', '*'],
    ['*', '*', '*', '0', '1', '*'],
    ['*', '*', '*', '*', '*', '*']
]

# Directions for checking (right, left, down, up, diagonals)
directions = [(0, 1), (0, -1), (1, 0), (-1, 0), (1, 1), (1, -1), (-1, 1), (-1, -1)]

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

# Round 1: Your move
board[4][3] = '0'
flip_pieces(board, 4, 3, '0')

# Round 2: Opponent's move
board[4][4] = '1'
flip_pieces(board, 4, 4, '1')

# Round 3: Your move
board[1][2] = '0'
flip_pieces(board, 1, 2, '0')

# Flatten the board for the final output
flattened_board = [cell for row in board for cell in row]
print(flattened_board)