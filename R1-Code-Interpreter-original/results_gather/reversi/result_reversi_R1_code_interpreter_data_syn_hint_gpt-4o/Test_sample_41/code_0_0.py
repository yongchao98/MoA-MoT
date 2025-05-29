# Initial board setup
board = [
    ['*', '0', '*', '*'],
    ['*', '0', '0', '*'],
    ['*', '0', '1', '0'],
    ['*', '*', '*', '1']
]

def flip_pieces(board, row, col, player):
    # Directions: right, left, down, up, diagonals
    directions = [(0, 1), (0, -1), (1, 0), (-1, 0), (1, 1), (1, -1), (-1, 1), (-1, -1)]
    opponent = '1' if player == '0' else '0'
    
    for dr, dc in directions:
        r, c = row + dr, col + dc
        pieces_to_flip = []
        
        while 0 <= r < 4 and 0 <= c < 4 and board[r][c] == opponent:
            pieces_to_flip.append((r, c))
            r += dr
            c += dc
        
        if 0 <= r < 4 and 0 <= c < 4 and board[r][c] == player:
            for rr, cc in pieces_to_flip:
                board[rr][cc] = player

# Round 1: Your move
board[2][3] = '0'
flip_pieces(board, 2, 3, '0')

# Round 2: Opponent's move
board[3][3] = '1'
flip_pieces(board, 3, 3, '1')

# Round 3: Your move
board[0][1] = '0'
flip_pieces(board, 0, 1, '0')

# Flatten the board to the required format
flattened_board = [cell for row in board for cell in row]
print(flattened_board)