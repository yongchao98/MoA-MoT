# Initial board setup
board = [
    ['0', '*', '*', '*'],
    ['*', '0', '1', '1'],
    ['*', '0', '0', '0'],
    ['*', '*', '*', '*']
]

# Function to place a piece and flip opponent's pieces
def place_piece(board, row, col, player):
    opponent = '1' if player == '0' else '0'
    board[row][col] = player
    
    # Directions: horizontal, vertical, and diagonal
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (-1, 1), (1, -1), (1, 1)]
    
    for dr, dc in directions:
        r, c = row + dr, col + dc
        pieces_to_flip = []
        
        while 0 <= r < 4 and 0 <= c < 4 and board[r][c] == opponent:
            pieces_to_flip.append((r, c))
            r += dr
            c += dc
        
        if 0 <= r < 4 and 0 <= c < 4 and board[r][c] == player:
            for rr, cc in pieces_to_flip:
                board[rr][cc] = player

# Round 1: Your move (Black) at (3,4)
place_piece(board, 2, 3, '0')

# Round 2: Opponent's move (White) at (2,4)
place_piece(board, 1, 3, '1')

# Round 3: Your move (Black) at (1,1)
place_piece(board, 0, 1, '0')

# Flatten the board to a single list
flattened_board = [cell for row in board for cell in row]
print(flattened_board)