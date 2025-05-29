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
        ['*', '*', '0', '*', '*', '*'],
        ['*', '*', '0', '0', '*', '*'],
        ['*', '1', '1', '0', '*', '*'],
        ['*', '*', '*', '*', '0', '*'],
        ['*', '*', '*', '*', '*', '*']
    ]
    
    # Round 1: Black places at (2,3)
    grid[1][2] = '0'
    flip_pieces(grid, 1, 2, '0')
    
    # Round 2: White places at (4,2)
    grid[3][1] = '1'
    flip_pieces(grid, 3, 1, '1')
    
    # Round 3: Black places at (5,5)
    grid[4][4] = '0'
    flip_pieces(grid, 4, 4, '0')
    
    # Flatten the grid for output
    result = [cell for row in grid for cell in row]
    formatted_result = ','.join(result)
    print(f"<<<{formatted_result}>>>")

play_game()