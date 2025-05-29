def flip_pieces(grid, x, y, player):
    # Directions: right, left, down, up, diagonals
    directions = [(0, 1), (0, -1), (1, 0), (-1, 0), (1, 1), (1, -1), (-1, 1), (-1, -1)]
    opponent = '1' if player == '0' else '0'
    
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        pieces_to_flip = []
        
        while 0 <= nx < 4 and 0 <= ny < 4 and grid[nx][ny] == opponent:
            pieces_to_flip.append((nx, ny))
            nx += dx
            ny += dy
        
        if 0 <= nx < 4 and 0 <= ny < 4 and grid[nx][ny] == player:
            for px, py in pieces_to_flip:
                grid[px][py] = player

# Initial grid setup
grid = [
    ['*', '*', '*', '*'],
    ['*', '1', '0', '*'],
    ['*', '0', '0', '*'],
    ['0', '0', '0', '*']
]

# Round 1: Your move (4,3)
grid[3][2] = '0'
flip_pieces(grid, 3, 2, '0')

# Round 2: Opponent's move (4,2)
grid[3][1] = '1'
flip_pieces(grid, 3, 1, '1')

# Round 3: Your move (4,1)
grid[3][0] = '0'
flip_pieces(grid, 3, 0, '0')

# Flatten the grid for the final output
flattened_grid = [cell for row in grid for cell in row]
print(flattened_grid)