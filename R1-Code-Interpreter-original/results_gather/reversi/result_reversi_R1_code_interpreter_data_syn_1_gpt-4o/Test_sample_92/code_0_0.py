def print_grid(grid):
    for row in grid:
        print(''.join(row))

def is_on_board(x, y, n):
    return 0 <= x < n and 0 <= y < n

def flip_pieces(grid, x, y, player, opponent, n):
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (-1, 1), (1, -1), (1, 1)]
    to_flip = []
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        pieces_to_flip = []
        while is_on_board(nx, ny, n) and grid[nx][ny] == opponent:
            pieces_to_flip.append((nx, ny))
            nx += dx
            ny += dy
        if is_on_board(nx, ny, n) and grid[nx][ny] == player:
            to_flip.extend(pieces_to_flip)
    for fx, fy in to_flip:
        grid[fx][fy] = player

def play_game(n, initial_grid, moves):
    grid = [list(row) for row in initial_grid]
    for i, (x, y) in enumerate(moves):
        player = '0' if i % 2 == 0 else '1'
        opponent = '1' if player == '0' else '0'
        grid[x][y] = player
        flip_pieces(grid, x, y, player, opponent, n)
    return grid

n = 6
initial_grid = [
    "******",
    "**0***",
    "**00**",
    "*110**",
    "****0*",
    "******"
]

moves = [(2, 3), (4, 2), (5, 5)]
final_grid = play_game(n, initial_grid, moves)
print_grid(final_grid)