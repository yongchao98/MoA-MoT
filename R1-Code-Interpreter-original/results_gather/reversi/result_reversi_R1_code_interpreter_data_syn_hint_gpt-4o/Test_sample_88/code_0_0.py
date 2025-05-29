def print_board(board):
    for row in board:
        print(''.join(row))

def is_on_board(x, y, n):
    return 0 <= x < n and 0 <= y < n

def flip_pieces(board, x, y, player, opponent, n):
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (-1, 1), (1, -1), (1, 1)]
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

def play_game(n, initial_board, moves):
    board = [list(row) for row in initial_board]
    player, opponent = '0', '1'
    
    for move in moves:
        x, y = move
        board[x][y] = player
        flip_pieces(board, x, y, player, opponent, n)
        player, opponent = opponent, player  # Switch turns
    
    return board

# Initial board configuration
initial_board = [
    "******",
    "******",
    "*000**",
    "**01**",
    "***01*",
    "******"
]

# Moves: (row, column) - 0-indexed
moves = [(4, 3), (4, 4), (2, 1)]

# Play the game
final_board = play_game(6, initial_board, moves)

# Print the final board
print_board(final_board)