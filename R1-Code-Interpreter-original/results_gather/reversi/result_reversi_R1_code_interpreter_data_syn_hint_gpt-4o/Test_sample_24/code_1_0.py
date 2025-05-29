def flip_pieces(grid, x, y, player):
    opponent = '1' if player == '0' else '0'
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
    grid = [
        ['*', '*', '*', '*'],
        ['*', '1', '1', '1'],
        ['*', '0', '0', '*'],
        ['*', '*', '0', '*']
    ]
    
    # Round 1: Black's move at (4,3)
    grid[3][2] = '0'
    flip_pieces(grid, 3, 2, '0')
    
    # Round 2: White's move at (2,4)
    grid[1][3] = '1'
    flip_pieces(grid, 1, 3, '1')
    
    # Flatten the grid for output
    result = [cell for row in grid for cell in row]
    print(result)

play_game()