def flip_pieces(board, row, col, player):
    opponent = '1' if player == '0' else '0'
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (-1, 1), (1, -1), (1, 1)]
    n = len(board)
    
    for dr, dc in directions:
        r, c = row + dr, col + dc
        pieces_to_flip = []
        
        while 0 <= r < n and 0 <= c < n and board[r][c] == opponent:
            pieces_to_flip.append((r, c))
            r += dr
            c += dc
        
        if 0 <= r < n and 0 <= c < n and board[r][c] == player:
            for rr, cc in pieces_to_flip:
                board[rr][cc] = player

def apply_move(board, row, col, player):
    board[row][col] = player
    flip_pieces(board, row, col, player)

def board_to_string(board):
    return ','.join([cell for row in board for cell in row])

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

# Output the final board state
print(board_to_string(board))