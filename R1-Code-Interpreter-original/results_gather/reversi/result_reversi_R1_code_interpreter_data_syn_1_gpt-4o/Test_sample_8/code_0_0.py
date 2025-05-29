def flip_pieces(board, x, y, player):
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (1, 1), (-1, 1), (1, -1)]
    n = len(board)
    opponent = '1' if player == '0' else '0'
    
    def is_valid(nx, ny):
        return 0 <= nx < n and 0 <= ny < n
    
    def flip_in_direction(dx, dy):
        nx, ny = x + dx, y + dy
        pieces_to_flip = []
        while is_valid(nx, ny) and board[nx][ny] == opponent:
            pieces_to_flip.append((nx, ny))
            nx += dx
            ny += dy
        if is_valid(nx, ny) and board[nx][ny] == player:
            for px, py in pieces_to_flip:
                board[px][py] = player
    
    for dx, dy in directions:
        flip_in_direction(dx, dy)

# Initial board setup
board = [
    ['*', '*', '*', '*'],
    ['*', '1', '0', '*'],
    ['*', '0', '0', '*'],
    ['*', '*', '0', '*']
]

# Place the black piece at (4, 3) which is (3, 2) in 0-indexed
board[3][2] = '0'
flip_pieces(board, 3, 2, '0')

# Flatten the board to the required output format
flattened_board = [cell for row in board for cell in row]
print(flattened_board)