def print_grid(grid):
    for row in grid:
        print(''.join(row))

def flip_pieces(grid, x, y, player):
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (-1, 1), (1, -1), (1, 1)]
    opponent = '1' if player == '0' else '0'
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
    grid = [
        ['*', '*', '*', '*', '*', '*'],
        ['*', '*', '*', '*', '0', '*'],
        ['*', '*', '1', '0', '1', '*'],
        ['*', '*', '0', '0', '*', '*'],
        ['*', '*', '*', '0', '*', '*'],
        ['*', '*', '*', '*', '*', '*']
    ]
    
    # Round 1: Your move (5,4) -> (4,3) in 0-based index
    grid[4][3] = '0'
    flip_pieces(grid, 4, 3, '0')
    
    # Round 2: Opponent's move (3,5) -> (2,4) in 0-based index
    grid[2][4] = '1'
    flip_pieces(grid, 2, 4, '1')
    
    # Round 3: Your move (2,5) -> (1,4) in 0-based index
    grid[1][4] = '0'
    flip_pieces(grid, 1, 4, '0')
    
    print_grid(grid)

play_game()