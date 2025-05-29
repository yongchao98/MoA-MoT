def flip_pieces(board, x, y, player):
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (-1, 1), (1, -1), (1, 1)]
    opponent = '1' if player == '0' else '0'
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
        ['*', '0', '*', '*'],
        ['*', '0', '0', '*'],
        ['*', '0', '1', '0'],
        ['*', '*', '*', '1']
    ]
    
    # Round 1: Your move (3,4)
    board[2][3] = '0'
    flip_pieces(board, 2, 3, '0')
    
    # Round 2: Opponent's move (4,4)
    board[3][3] = '1'
    flip_pieces(board, 3, 3, '1')
    
    # Round 3: Your move (1,2)
    board[0][1] = '0'
    flip_pieces(board, 0, 1, '0')
    
    # Flatten the board to a single list
    flat_board = [cell for row in board for cell in row]
    print(flat_board)

play_game()