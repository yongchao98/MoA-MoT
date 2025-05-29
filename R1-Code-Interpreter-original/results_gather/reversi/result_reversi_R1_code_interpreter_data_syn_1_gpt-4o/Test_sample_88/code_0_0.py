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
        ['*', '*', '*', '*', '*', '*'],
        ['*', '*', '*', '*', '*', '*'],
        ['*', '0', '0', '0', '*', '*'],
        ['*', '*', '0', '1', '*', '*'],
        ['*', '*', '*', '0', '1', '*'],
        ['*', '*', '*', '*', '*', '*']
    ]
    
    # Round 1: Your move (Black)
    grid[4][3] = '0'
    flip_pieces(grid, '0', '1', 4, 3)
    
    # Round 2: Opponent's move (White)
    grid[4][4] = '1'
    flip_pieces(grid, '1', '0', 4, 4)
    
    # Round 3: Your move (Black)
    grid[2][1] = '0'
    flip_pieces(grid, '0', '1', 2, 1)
    
    return grid

final_grid = play_game()
print_grid(final_grid)