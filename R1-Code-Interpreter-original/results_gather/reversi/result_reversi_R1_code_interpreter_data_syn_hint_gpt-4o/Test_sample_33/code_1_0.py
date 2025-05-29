# Initial board setup
board = [
    ['*', '*', '*', '*'],
    ['*', '1', '0', '*'],
    ['*', '0', '1', '0'],
    ['*', '*', '*', '1']
]

# Function to flip pieces
def flip_pieces(board, x, y, player):
    opponent = '1' if player == '0' else '0'
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (-1, 1), (1, -1), (1, 1)]
    
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

# Round 1: Your move
board[2][3] = '0'  # Place black piece at (3,4)
flip_pieces(board, 2, 3, '0')

# Round 2: Opponent's move
board[3][3] = '1'  # Place white piece at (4,4)
flip_pieces(board, 3, 3, '1')

# Flatten the board to the required format
flattened_board = [cell for row in board for cell in row]
print(flattened_board)