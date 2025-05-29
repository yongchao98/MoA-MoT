def flip_pieces(board, row, col, player):
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (-1, 1), (1, -1), (1, 1)]
    opponent = '1' if player == '0' else '0'
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

def apply_move(board, move, player):
    row, col = move
    board[row][col] = player
    flip_pieces(board, row, col, player)

# Initial board setup
board = [
    ['*', '*', '*', '*', '*', '*'],
    ['*', '*', '*', '1', '*', '*'],
    ['*', '0', '0', '1', '*', '*'],
    ['*', '*', '0', '0', '0', '*'],
    ['*', '*', '*', '*', '*', '*'],
    ['*', '*', '*', '*', '*', '*']
]

# Apply moves
apply_move(board, (3, 2), '0')  # Round 1: Your move
apply_move(board, (2, 4), '1')  # Round 2: Opponent's move
apply_move(board, (4, 5), '0')  # Round 3: Your move

# Flatten the board for the final output
flattened_board = [cell for row in board for cell in row]
print(flattened_board)