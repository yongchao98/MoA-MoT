def print_board(board):
    return ','.join([','.join(row) for row in board])

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

def play_game():
    board = [
        ['*', '*', '*', '*', '*', '*'],
        ['*', '*', '*', '*', '*', '*'],
        ['*', '*', '1', '0', '*', '*'],
        ['*', '*', '0', '1', '*', '*'],
        ['*', '*', '*', '0', '1', '*'],
        ['*', '*', '*', '*', '*', '*']
    ]
    
    # Round 1: Black plays at (5, 4)
    board[4][3] = '0'
    flip_pieces(board, 4, 3, '0')
    
    # Round 2: White plays at (5, 5)
    board[4][4] = '1'
    flip_pieces(board, 4, 4, '1')
    
    return print_board(board)

print(play_game())