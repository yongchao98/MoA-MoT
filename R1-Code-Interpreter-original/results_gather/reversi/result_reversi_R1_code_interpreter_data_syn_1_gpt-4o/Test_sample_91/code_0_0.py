def print_board(board):
    for row in board:
        print(''.join(row))

def flip_pieces(board, x, y, player):
    opponent = '1' if player == '0' else '0'
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (-1, 1), (1, -1), (1, 1)]
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

def apply_move(board, x, y, player):
    board[x][y] = player
    flip_pieces(board, x, y, player)

# Initial board configuration
board = [
    ['*', '*', '*', '*', '*', '*'],
    ['*', '1', '*', '*', '*', '*'],
    ['*', '0', '1', '0', '*', '*'],
    ['*', '*', '0', '0', '*', '*'],
    ['*', '*', '*', '0', '*', '*'],
    ['*', '*', '*', '*', '*', '*']
]

# Apply moves
apply_move(board, 3, 2, '0')  # Round 1
apply_move(board, 2, 2, '1')  # Round 2
apply_move(board, 5, 4, '0')  # Round 3

# Print final board
print_board(board)