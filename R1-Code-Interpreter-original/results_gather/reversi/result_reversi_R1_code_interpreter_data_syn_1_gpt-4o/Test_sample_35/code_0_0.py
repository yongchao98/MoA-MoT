def flip_pieces(board, row, col, player):
    n = len(board)
    opponent = '1' if player == '0' else '0'
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (1, 1), (-1, 1), (1, -1)]
    
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

# Initial board state
board = [
    ['*', '*', '*', '*'],
    ['*', '1', '0', '*'],
    ['*', '0', '1', '*'],
    ['*', '*', '0', '1']
]

# Apply Round 1: Black plays at (4,3)
apply_move(board, 3, 2, '0')

# Apply Round 2: White plays at (4,4)
apply_move(board, 3, 3, '1')

# Flatten the board to a single list
final_board = [cell for row in board for cell in row]

print(final_board)