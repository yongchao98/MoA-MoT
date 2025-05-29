def print_grid(grid):
    for row in grid:
        print(''.join(row))

def is_valid_position(x, y, n):
    return 0 <= x < n and 0 <= y < n

def flip_pieces(grid, x, y, player, opponent):
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1), 
                  (-1, -1), (-1, 1), (1, -1), (1, 1)]
    n = len(grid)
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        pieces_to_flip = []
        while is_valid_position(nx, ny, n) and grid[nx][ny] == opponent:
            pieces_to_flip.append((nx, ny))
            nx += dx
            ny += dy
        if is_valid_position(nx, ny, n) and grid[nx][ny] == player:
            for px, py in pieces_to_flip:
                grid[px][py] = player

def play_game():
    n = 6
    grid = [
        ['*', '*', '*', '*', '*', '*'],
        ['*', '*', '*', '*', '*', '*'],
        ['*', '0', '0', '0', '*', '*'],
        ['*', '0', '1', '1', '*', '*'],
        ['0', '*', '*', '*', '*', '*'],
        ['*', '*', '*', '*', '*', '*']
    ]

    # Round 1: Your's landing position is (3,2)
    grid[2][1] = '0'
    flip_pieces(grid, 2, 1, '0', '1')

    # Round 2: The opponent's landing position is (4,2)
    grid[3][1] = '1'
    flip_pieces(grid, 3, 1, '1', '0')

    # Round 3: Your's landing position is (5,1)
    grid[4][0] = '0'
    flip_pieces(grid, 4, 0, '0', '1')

    print_grid(grid)

play_game()