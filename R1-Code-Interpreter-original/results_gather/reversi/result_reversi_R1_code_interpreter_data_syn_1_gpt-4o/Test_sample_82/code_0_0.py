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

def play_game():
    board = [
        ['*', '*', '*', '*', '*', '*'],
        ['*', '*', '*', '*', '*', '*'],
        ['*', '0', '0', '0', '*', '*'],
        ['*', '0', '1', '1', '*', '*'],
        ['0', '*', '*', '*', '*', '*'],
        ['*', '*', '*', '*', '*', '*']
    ]
    
    # Round 1: Your move (Black) at (3,2)
    board[2][1] = '0'
    flip_pieces(board, 2, 1, '0')
    
    # Round 2: Opponent's move (White) at (4,2)
    board[3][1] = '1'
    flip_pieces(board, 3, 1, '1')
    
    # Round 3: Your move (Black) at (5,1)
    board[4][0] = '0'
    flip_pieces(board, 4, 0, '0')
    
    print_board(board)

play_game()