def print_board(board):
    for row in board:
        print(''.join(row))

def flip_pieces(board, player, opponent, x, y):
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
        ['*', '*', '0', '*'],
        ['*', '1', '0', '1'],
        ['*', '0', '0', '*'],
        ['*', '*', '0', '*']
    ]
    
    # Round 1: Your move (Black) at (4,3)
    board[3][2] = '0'
    flip_pieces(board, '0', '1', 3, 2)
    
    # Round 2: Opponent's move (White) at (2,4)
    board[1][3] = '1'
    flip_pieces(board, '1', '0', 1, 3)
    
    # Round 3: Your move (Black) at (1,3)
    board[0][2] = '0'
    flip_pieces(board, '0', '1', 0, 2)
    
    # Print the final board
    print_board(board)

play_game()