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
        ['*', '*', '*', '*'],
        ['0', '0', '0', '*'],
        ['1', '1', '1', '*'],
        ['*', '*', '*', '*']
    ]
    
    # Round 1: Your move (Black) at (2, 1)
    board[1][0] = '0'
    flip_pieces(board, 1, 0, '0')
    
    # Round 2: Opponent's move (White) at (3, 1)
    board[2][0] = '1'
    flip_pieces(board, 2, 0, '1')
    
    # Flatten the board to a single list
    flat_board = [cell for row in board for cell in row]
    return flat_board

result = play_game()
print(result)