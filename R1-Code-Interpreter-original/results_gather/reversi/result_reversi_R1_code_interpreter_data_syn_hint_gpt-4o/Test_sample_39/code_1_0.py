def flip_pieces(board, row, col, player):
    # Directions: horizontal, vertical, and diagonals
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1), 
                  (-1, -1), (-1, 1), (1, -1), (1, 1)]
    opponent = '1' if player == '0' else '0'
    n = len(board)
    
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
    # Initial board setup
    board = [
        ['*', '*', '*', '*'],
        ['*', '1', '1', '1'],
        ['*', '0', '0', '0'],
        ['*', '*', '*', '*']
    ]
    
    # Round 1: Black's move at (3,4)
    board[2][3] = '0'
    flip_pieces(board, 2, 3, '0')
    
    # Round 2: White's move at (2,4)
    board[1][3] = '1'
    flip_pieces(board, 1, 3, '1')
    
    # Flatten the board to a single list
    flat_board = [cell for row in board for cell in row]
    return flat_board

# Execute the game simulation
result = play_game()
print(result)