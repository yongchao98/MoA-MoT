def print_board(board):
    for row in board:
        print(''.join(row))

def is_valid_position(x, y, n):
    return 0 <= x < n and 0 <= y < n

def flip_pieces(board, x, y, player):
    n = len(board)
    opponent = '1' if player == '0' else '0'
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (-1, 1), (1, -1), (1, 1)]
    
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        pieces_to_flip = []
        
        while is_valid_position(nx, ny, n) and board[nx][ny] == opponent:
            pieces_to_flip.append((nx, ny))
            nx += dx
            ny += dy
        
        if is_valid_position(nx, ny, n) and board[nx][ny] == player:
            for px, py in pieces_to_flip:
                board[px][py] = player

def play_game():
    n = 6
    board = [
        ['*', '*', '*', '*', '*', '*'],
        ['0', '0', '0', '*', '*', '*'],
        ['*', '*', '1', '0', '*', '*'],
        ['*', '*', '0', '1', '*', '*'],
        ['*', '*', '*', '*', '*', '*'],
        ['*', '*', '*', '*', '*', '*']
    ]
    
    # Round 1: Your move (Black) at (2, 3)
    board[1][2] = '0'
    flip_pieces(board, 1, 2, '0')
    
    # Round 2: Opponent's move (White) at (2, 2)
    board[1][1] = '1'
    flip_pieces(board, 1, 1, '1')
    
    # Round 3: Your move (Black) at (2, 1)
    board[1][0] = '0'
    flip_pieces(board, 1, 0, '0')
    
    print_board(board)

play_game()