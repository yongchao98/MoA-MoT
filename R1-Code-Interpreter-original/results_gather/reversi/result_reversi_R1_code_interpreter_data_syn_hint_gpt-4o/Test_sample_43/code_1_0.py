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
        ['*', '0', '*', '*'],
        ['*', '0', '0', '*'],
        ['1', '1', '0', '*'],
        ['*', '*', '*', '0']
    ]
    
    # Round 1: Black plays at (1,2)
    grid[0][1] = '0'
    flip_pieces(grid, 0, 1, '0')
    
    # Round 2: White plays at (3,1)
    grid[2][0] = '1'
    flip_pieces(grid, 2, 0, '1')
    
    # Round 3: Black plays at (4,4)
    grid[3][3] = '0'
    flip_pieces(grid, 3, 3, '0')
    
    # Flatten the grid for output
    flat_grid = [cell for row in grid for cell in row]
    print(flat_grid)

play_game()