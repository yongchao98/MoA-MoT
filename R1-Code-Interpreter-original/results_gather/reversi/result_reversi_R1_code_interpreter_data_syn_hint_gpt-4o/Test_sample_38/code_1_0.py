def flip_pieces(board, x, y, player):
    # Directions: right, left, down, up, down-right, down-left, up-right, up-left
    directions = [(0, 1), (0, -1), (1, 0), (-1, 0), (1, 1), (1, -1), (-1, 1), (-1, -1)]
    opponent = '1' if player == '0' else '0'
    
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        pieces_to_flip = []
        
        while 0 <= nx < 4 and 0 <= ny < 4 and board[nx][ny] == opponent:
            pieces_to_flip.append((nx, ny))
            nx += dx
            ny += dy
        
        if 0 <= nx < 4 and 0 <= ny < 4 and board[nx][ny] == player:
            for px, py in pieces_to_flip:
                board[px][py] = player

# Initial board setup
board = [
    ['*', '*', '*', '*'],
    ['*', '1', '1', '1'],
    ['*', '0', '0', '*'],
    ['*', '*', '0', '*']
]

# Round 1: Black's move at (4,3) -> (3,2) in 0-indexed
board[3][2] = '0'
flip_pieces(board, 3, 2, '0')

# Round 2: White's move at (2,4) -> (1,3) in 0-indexed
board[1][3] = '1'
flip_pieces(board, 1, 3, '1')

# Flatten the board to the required output format
flattened_board = [cell for row in board for cell in row]
print(flattened_board)