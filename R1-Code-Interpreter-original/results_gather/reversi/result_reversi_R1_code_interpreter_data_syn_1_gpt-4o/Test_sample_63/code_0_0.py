# Initial board setup
board = [
    ['*', '*', '*', '*', '*', '*'],
    ['*', '*', '0', '*', '*', '*'],
    ['*', '*', '0', '0', '*', '*'],
    ['*', '1', '1', '1', '*', '*'],
    ['*', '*', '*', '*', '*', '*'],
    ['*', '*', '*', '*', '*', '*']
]

# Function to flip pieces
def flip_pieces(board, x, y, player):
    opponent = '1' if player == '0' else '0'
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (1, 1), (-1, 1), (1, -1)]
    n = len(board)
    
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

# Round 1: Black places at (2,3)
board[1][2] = '0'
flip_pieces(board, 1, 2, '0')

# Round 2: White places at (4,2)
board[3][1] = '1'
flip_pieces(board, 3, 1, '1')

# Convert board to the required format
result = [cell for row in board for cell in row]
print(result)