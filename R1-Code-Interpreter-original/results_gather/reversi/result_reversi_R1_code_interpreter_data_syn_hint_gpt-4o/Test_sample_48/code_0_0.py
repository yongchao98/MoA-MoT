def print_grid(grid):
    for row in grid:
        print(''.join(row))

def is_valid_position(x, y, n):
    return 0 <= x < n and 0 <= y < n

def flip_pieces(grid, x, y, player):
    n = len(grid)
    opponent = '1' if player == '0' else '0'
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (-1, 1), (1, -1), (1, 1)]
    
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
    n = 4
    grid = [
        ['*', '*', '*', '*'],
        ['0', '0', '0', '*'],
        ['1', '1', '0', '*'],
        ['*', '*', '*', '0']
    ]
    
    # Round 1: Your move (Black) at (2, 1)
    grid[1][0] = '0'
    flip_pieces(grid, 1, 0, '0')
    
    # Round 2: Opponent's move (White) at (3, 1)
    grid[2][0] = '1'
    flip_pieces(grid, 2, 0, '1')
    
    # Round 3: Your move (Black) at (4, 4)
    grid[3][3] = '0'
    flip_pieces(grid, 3, 3, '0')
    
    # Print the final grid
    print_grid(grid)

play_game()