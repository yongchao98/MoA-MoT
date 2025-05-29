# Initial board setup
board = [
    ['*', '*', '*', '*'],
    ['*', '1', '1', '1'],
    ['*', '0', '0', '*'],
    ['*', '*', '0', '*']
]

# Directions for checking: (row_offset, col_offset)
directions = [
    (0, 1), (1, 0), (0, -1), (-1, 0),  # Horizontal and vertical
    (1, 1), (1, -1), (-1, 1), (-1, -1) # Diagonals
]

def is_valid_position(x, y, n):
    return 0 <= x < n and 0 <= y < n

def flip_pieces(board, x, y, player):
    opponent = '0' if player == '1' else '1'
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

# Round 1: Black's move at (4,3) -> (3,2) in 0-indexed
board[3][2] = '0'
flip_pieces(board, 3, 2, '0')

# Round 2: White's move at (2,4) -> (1,3) in 0-indexed
board[1][3] = '1'
flip_pieces(board, 1, 3, '1')

# Flatten the board for the final output
flattened_board = [cell for row in board for cell in row]
print(flattened_board)