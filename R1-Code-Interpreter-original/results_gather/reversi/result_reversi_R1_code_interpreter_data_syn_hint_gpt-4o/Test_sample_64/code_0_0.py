# Initial board configuration
board = [
    ['*', '*', '*', '*', '*', '*'],
    ['*', '*', '0', '1', '*', '*'],
    ['*', '*', '0', '1', '*', '*'],
    ['*', '*', '0', '1', '*', '*'],
    ['*', '*', '*', '*', '*', '*'],
    ['*', '*', '*', '*', '*', '*']
]

# Function to flip pieces
def flip_pieces(board, x, y, player):
    opponent = '1' if player == '0' else '0'
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (1, 1), (-1, 1), (1, -1)]
    
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

# Round 1: Your move (black) at (2, 3)
board[1][2] = '0'
flip_pieces(board, 1, 2, '0')

# Round 2: Opponent's move (white) at (2, 4)
board[1][3] = '1'
flip_pieces(board, 1, 3, '1')

# Print the final board configuration
final_board = [cell for row in board for cell in row]
print(final_board)