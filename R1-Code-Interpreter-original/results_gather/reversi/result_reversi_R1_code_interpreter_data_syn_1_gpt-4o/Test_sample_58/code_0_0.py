def apply_move(board, move, player):
    n = len(board)
    x, y = move
    board[x][y] = player
    
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (-1, 1), (1, -1), (1, 1)]
    
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        pieces_to_flip = []
        
        while 0 <= nx < n and 0 <= ny < n:
            if board[nx][ny] == '*':
                break
            if board[nx][ny] == player:
                for px, py in pieces_to_flip:
                    board[px][py] = player
                break
            pieces_to_flip.append((nx, ny))
            nx += dx
            ny += dy

# Initial board
board = [
    ['*', '0', '*', '*'],
    ['*', '0', '0', '*'],
    ['*', '0', '1', '0'],
    ['*', '*', '*', '1']
]

# Apply moves
apply_move(board, (2, 3), '0')  # Round 1: Your move
apply_move(board, (3, 3), '1')  # Round 2: Opponent's move
apply_move(board, (0, 1), '0')  # Round 3: Your move

# Flatten the board for the final output
flattened_board = [cell for row in board for cell in row]
print(flattened_board)