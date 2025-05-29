def flip_pieces(board, row, col, player):
    opponent = '1' if player == '0' else '0'
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1), 
                  (-1, -1), (-1, 1), (1, -1), (1, 1)]
    
    for dr, dc in directions:
        r, c = row + dr, col + dc
        pieces_to_flip = []
        while 0 <= r < 6 and 0 <= c < 6 and board[r][c] == opponent:
            pieces_to_flip.append((r, c))
            r += dr
            c += dc
        if 0 <= r < 6 and 0 <= c < 6 and board[r][c] == player:
            for rr, cc in pieces_to_flip:
                board[rr][cc] = player

def print_board(board):
    for row in board:
        print(''.join(row))

# Initial board setup
board = [
    ['*', '*', '*', '*', '*', '*'],
    ['*', '*', '*', '*', '*', '*'],
    ['*', '*', '1', '1', '1', '*'],
    ['*', '*', '0', '0', '*', '*'],
    ['*', '*', '*', '0', '*', '*'],
    ['*', '*', '*', '*', '*', '*']
]

# Round 1: Black's move
board[4][3] = '0'
flip_pieces(board, 4, 3, '0')

# Round 2: White's move
board[2][4] = '1'
flip_pieces(board, 2, 4, '1')

# Print final board state
print_board(board)