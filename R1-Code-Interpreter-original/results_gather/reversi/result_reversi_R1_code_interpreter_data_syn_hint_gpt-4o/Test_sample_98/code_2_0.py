def print_grid(grid):
    for row in grid:
        print(''.join(row))

def flip_pieces(grid, player, opponent, x, y):
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (1, 1), (-1, 1), (1, -1)]
    n = len(grid)
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        pieces_to_flip = []
        while 0 <= nx < n and 0 <= ny < n and grid[nx][ny] == opponent:
            pieces_to_flip.append((nx, ny))
            nx += dx
            ny += dy
        if 0 <= nx < n and 0 <= ny < n and grid[nx][ny] == player:
            for px, py in pieces_to_flip:
                grid[px][py] = player

def play_game():
    n = 6
    grid = [
        ['*', '0', '*', '*', '*', '*'],
        ['*', '0', '*', '*', '*', '*'],
        ['*', '0', '1', '0', '*', '*'],
        ['*', '*', '0', '1', '*', '*'],
        ['*', '*', '*', '*', '*', '*'],
        ['*', '*', '*', '*', '*', '*']
    ]

    # Round 1: Your move (Black) at (3, 2)
    grid[2][1] = '0'
    flip_pieces(grid, '0', '1', 2, 1)

    # Round 2: Opponent's move (White) at (2, 2)
    grid[1][1] = '1'
    flip_pieces(grid, '1', '0', 1, 1)

    # Round 3: Your move (Black) at (1, 2)
    grid[0][1] = '0'
    flip_pieces(grid, '0', '1', 0, 1)

    # Print the final grid
    print_grid(grid)

play_game()