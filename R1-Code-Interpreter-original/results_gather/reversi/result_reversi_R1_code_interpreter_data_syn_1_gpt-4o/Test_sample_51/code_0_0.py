def flip_pieces(board, row, col, player):
    n = len(board)
    opponent = '1' if player == '0' else '0'
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (1, 1), (-1, 1), (1, -1)]
    
    for dr, dc in directions:
        r, c = row + dr, col + dc
        pieces_to_flip = []
        
        while 0 <= r < n and 0 <= c < n and board[r][c] == opponent:
            pieces_to_flip.append((r, c))
            r += dr
            c += dc
        
        if 0 <= r < n and 0 <= c < n and board[r][c] == player:
            for rr, cc in pieces_to_flip:
                board[rr][cc] = player

def play_game():
    board = [
        ['*', '*', '*', '0'],
        ['*', '1', '0', '1'],
        ['*', '0', '0', '*'],
        ['*', '*', '0', '*']
    ]
    
    # Round 1: Your move (Black) at (4,3)
    board[3][2] = '0'
    flip_pieces(board, 3, 2, '0')
    
    # Round 2: Opponent's move (White) at (2,4)
    board[1][3] = '1'
    flip_pieces(board, 1, 3, '1')
    
    # Round 3: Your move (Black) at (1,4)
    board[0][3] = '0'
    flip_pieces(board, 0, 3, '0')
    
    # Flatten the board to a single list
    flat_board = [cell for row in board for cell in row]
    return flat_board

result = play_game()
print(result)