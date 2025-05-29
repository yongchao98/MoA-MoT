def print_board(board):
    for row in board:
        print(''.join(row))

def flip_pieces(board, x, y, player):
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (-1, 1), (1, -1), (1, 1)]
    opponent = '1' if player == '0' else '0'
    n = len(board)
    
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        pieces_to_flip = []
        
        while 0 <= nx < n and 0 <= ny < n:
            if board[nx][ny] == opponent:
                pieces_to_flip.append((nx, ny))
            elif board[nx][ny] == player:
                for px, py in pieces_to_flip:
                    board[px][py] = player
                break
            else:
                break
            nx += dx
            ny += dy

# Initial board setup
board = [
    ['1', '*', '*', '*'],
    ['0', '1', '0', '*'],
    ['*', '0', '1', '*'],
    ['*', '*', '*', '*']
]

# Round 1: Your move (Black) at (2,1)
board[1][0] = '0'
flip_pieces(board, 1, 0, '0')

# Round 2: Opponent's move (White) at (1,1)
board[0][0] = '1'
flip_pieces(board, 0, 0, '1')

# Print the final board state
print_board(board)