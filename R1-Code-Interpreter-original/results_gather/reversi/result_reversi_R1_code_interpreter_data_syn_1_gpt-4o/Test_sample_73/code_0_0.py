def print_board(board):
    for row in board:
        print(''.join(row))

def flip_pieces(board, x, y, player):
    opponent = '1' if player == '0' else '0'
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1), 
                  (-1, -1), (-1, 1), (1, -1), (1, 1)]
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
    n = 6
    board = [['*'] * n for _ in range(n)]
    board[2][2] = '0'
    board[2][3] = '1'
    board[3][2] = '1'
    board[3][3] = '0'
    board[4][2] = '1'
    
    # Round 1: Black's move
    board[2][3] = '0'
    flip_pieces(board, 2, 3, '0')
    
    # Round 2: White's move
    board[2][4] = '1'
    flip_pieces(board, 2, 4, '1')
    
    print_board(board)

play_game()