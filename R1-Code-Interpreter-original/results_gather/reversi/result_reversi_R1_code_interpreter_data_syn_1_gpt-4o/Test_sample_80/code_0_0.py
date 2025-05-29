def print_board(board):
    for row in board:
        print(''.join(row))

def apply_move(board, row, col, player):
    opponent = '1' if player == '0' else '0'
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (-1, 1), (1, -1), (1, 1)]
    board[row][col] = player
    
    for dr, dc in directions:
        r, c = row + dr, col + dc
        pieces_to_flip = []
        while 0 <= r < len(board) and 0 <= c < len(board) and board[r][c] == opponent:
            pieces_to_flip.append((r, c))
            r += dr
            c += dc
        if 0 <= r < len(board) and 0 <= c < len(board) and board[r][c] == player:
            for rr, cc in pieces_to_flip:
                board[rr][cc] = player

# Initial board setup
board = [
    ['*', '*', '*', '*', '*', '*'],
    ['*', '1', '*', '*', '*', '*'],
    ['*', '0', '1', '0', '*', '*'],
    ['*', '*', '0', '0', '0', '*'],
    ['*', '*', '*', '*', '*', '*'],
    ['*', '*', '*', '*', '*', '*']
]

# Apply moves
apply_move(board, 3, 2, '0')  # Round 1
apply_move(board, 2, 2, '1')  # Round 2
apply_move(board, 4, 5, '0')  # Round 3

# Print final board
print_board(board)