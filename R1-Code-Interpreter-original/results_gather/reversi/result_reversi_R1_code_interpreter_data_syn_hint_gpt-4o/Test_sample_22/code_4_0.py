# Initial board setup
board = [
    ['*', '*', '*', '*'],
    ['*', '1', '1', '1'],
    ['*', '0', '0', '0'],
    ['*', '*', '*', '*']
]

# Directions for checking: (row_change, col_change)
directions = [
    (0, 1), (0, -1), (1, 0), (-1, 0),  # Horizontal and vertical
    (1, 1), (1, -1), (-1, 1), (-1, -1) # Diagonals
]

def is_valid_position(x, y, n):
    return 0 <= x < n and 0 <= y < n

def flip_pieces(board, x, y, player):
    opponent = '1' if player == '0' else '0'
    n = len(board)
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        pieces_to_flip = []
        while is_valid_position(nx, ny, n) and board[nx][ny] == opponent:
            pieces_to_flip.append((nx, ny))
            nx += dx
            ny += dy
        if is_valid_position(nx, ny, n) and board[nx][ny] == player:
            for px, py in pieces_to_flip:
                board[px][py] = player

# Round 1: Black's move at (3, 4)
board[2][3] = '0'
flip_pieces(board, 2, 3, '0')

# Round 2: White's move at (2, 4)
board[1][3] = '1'
flip_pieces(board, 1, 3, '1')

# Flatten the board to the required output format
flattened_board = [cell for row in board for cell in row]
print(flattened_board)