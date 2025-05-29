def print_board(board):
    for row in board:
        print(''.join(row))

def is_on_board(x, y, n):
    return 0 <= x < n and 0 <= y < n

def flip_pieces(board, x, y, player, opponent, n):
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1), 
                  (-1, -1), (-1, 1), (1, -1), (1, 1)]
    to_flip = []
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        pieces_to_flip = []
        while is_on_board(nx, ny, n) and board[nx][ny] == opponent:
            pieces_to_flip.append((nx, ny))
            nx += dx
            ny += dy
        if is_on_board(nx, ny, n) and board[nx][ny] == player:
            to_flip.extend(pieces_to_flip)
    for fx, fy in to_flip:
        board[fx][fy] = player

def play_game():
    n = 4
    board = [
        ['*', '0', '*', '*'],
        ['*', '0', '0', '*'],
        ['*', '0', '1', '*'],
        ['*', '*', '0', '1']
    ]
    
    # Round 1: Your move (Black) at (4,3)
    board[3][2] = '0'
    flip_pieces(board, 3, 2, '0', '1', n)
    
    # Round 2: Opponent's move (White) at (4,4)
    board[3][3] = '1'
    flip_pieces(board, 3, 3, '1', '0', n)
    
    # Round 3: Your move (Black) at (1,2)
    board[0][1] = '0'
    flip_pieces(board, 0, 1, '0', '1', n)
    
    print_board(board)

play_game()