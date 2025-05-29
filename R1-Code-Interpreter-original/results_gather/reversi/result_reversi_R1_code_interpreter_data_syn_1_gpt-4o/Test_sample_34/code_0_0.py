def flip_pieces(board, x, y, player):
    n = len(board)
    opponent = '1' if player == '0' else '0'
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (1, 1), (-1, 1), (1, -1)]
    
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        pieces_to_flip = []
        
        while 0 <= nx < n and 0 <= ny < n and board[nx][ny] == opponent:
            pieces_to_flip.append((nx, ny))
            nx += dx
            ny += dy
        
        if 0 <= nx < n and 0 <= ny < n and board[nx][ny] == player:
            for px, py in pieces_to_flip:
                board[px][py] = player

def apply_move(board, x, y, player):
    board[x][y] = player
    flip_pieces(board, x, y, player)

# Initial board state
board = [
    ['*', '0', '*', '*'],
    ['*', '0', '0', '*'],
    ['1', '1', '1', '*'],
    ['*', '*', '*', '*']
]

# Round 1: Black's move
apply_move(board, 0, 1, '0')

# Round 2: White's move
apply_move(board, 2, 0, '1')

# Convert board to the required format
result = [cell for row in board for cell in row]
print(result)