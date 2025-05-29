def print_board(board):
    for row in board:
        print(''.join(row))

def is_on_board(x, y, n):
    return 0 <= x < n and 0 <= y < n

def flip_pieces(board, x, y, player, opponent, n):
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (-1, 1), (1, -1), (1, 1)]
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        pieces_to_flip = []
        while is_on_board(nx, ny, n) and board[nx][ny] == opponent:
            pieces_to_flip.append((nx, ny))
            nx += dx
            ny += dy
        if is_on_board(nx, ny, n) and board[nx][ny] == player:
            for fx, fy in pieces_to_flip:
                board[fx][fy] = player

def play_game(n, initial_board, moves):
    board = [list(row) for row in initial_board]
    for i, (x, y) in enumerate(moves):
        player = '0' if i % 2 == 0 else '1'
        opponent = '1' if player == '0' else '0'
        board[x][y] = player
        flip_pieces(board, x, y, player, opponent, n)
    return board

n = 6
initial_board = [
    "******",
    "**0***",
    "**00**",
    "**01**",
    "***01*",
    "******"
]

moves = [(4, 3), (4, 4), (1, 2)]  # 0-indexed positions

final_board = play_game(n, initial_board, moves)
print_board(final_board)