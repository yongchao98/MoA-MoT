# Initial board setup
board = [
    ['*', '*', '*', '*', '*', '*'],
    ['*', '1', '*', '*', '*', '*'],
    ['*', '0', '1', '0', '*', '*'],
    ['*', '*', '0', '1', '*', '*'],
    ['*', '*', '*', '*', '*', '*'],
    ['*', '*', '*', '*', '*', '*']
]

# Function to flip pieces
def flip_pieces(board, x, y, player):
    opponent = '1' if player == '0' else '0'
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (-1, 1), (1, -1), (1, 1)]
    
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        pieces_to_flip = []
        
        while 0 <= nx < 6 and 0 <= ny < 6 and board[nx][ny] == opponent:
            pieces_to_flip.append((nx, ny))
            nx += dx
            ny += dy
        
        if 0 <= nx < 6 and 0 <= ny < 6 and board[nx][ny] == player:
            for px, py in pieces_to_flip:
                board[px][py] = player

# Round 1: Your move (black) at (3, 2)
board[2][1] = '0'
flip_pieces(board, 2, 1, '0')

# Round 2: Opponent's move (white) at (2, 2)
board[1][1] = '1'
flip_pieces(board, 1, 1, '1')

# Flatten the board to the required format
flattened_board = [cell for row in board for cell in row]
print(flattened_board)