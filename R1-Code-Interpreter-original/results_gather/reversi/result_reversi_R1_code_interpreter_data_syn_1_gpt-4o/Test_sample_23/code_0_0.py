# Initial board setup
board = [
    ['1', '*', '*', '*'],
    ['0', '1', '0', '*'],
    ['*', '0', '1', '*'],
    ['*', '*', '*', '*']
]

# Directions for checking: (row_change, col_change)
directions = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (-1, 1), (1, -1), (1, 1)]

def flip_pieces(board, row, col, player):
    opponent = '0' if player == '1' else '1'
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

# Round 1: Black's move at (2, 1)
board[1][0] = '0'
flip_pieces(board, 1, 0, '0')

# Round 2: White's move at (1, 1)
board[0][0] = '1'
flip_pieces(board, 0, 0, '1')

# Flatten the board to a single list
flattened_board = [cell for row in board for cell in row]
print(flattened_board)