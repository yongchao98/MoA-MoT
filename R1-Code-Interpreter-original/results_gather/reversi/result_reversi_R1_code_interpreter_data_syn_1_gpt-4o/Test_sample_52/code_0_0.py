def print_board(board):
    for row in board:
        print(' '.join(row))
    print()

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
        ['1', '0', '*', '0'],
        ['0', '0', '0', '*'],
        ['*', '0', '1', '*'],
        ['*', '*', '*', '*']
    ]
    
    # Round 1: Your move (Black) at (2,1)
    board[1][0] = '0'
    flip_pieces(board, 1, 0, '0')
    
    # Round 2: Opponent's move (White) at (1,3)
    board[0][2] = '1'
    flip_pieces(board, 0, 2, '1')
    
    # Round 3: Your move (Black) at (1,4)
    board[0][3] = '0'
    flip_pieces(board, 0, 3, '0')
    
    # Flatten the board to a single list
    flat_board = [cell for row in board for cell in row]
    print(flat_board)

play_game()