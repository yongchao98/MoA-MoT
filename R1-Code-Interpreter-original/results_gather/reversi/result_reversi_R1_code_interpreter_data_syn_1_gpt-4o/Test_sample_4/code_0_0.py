def flip_pieces(grid, x, y, player):
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1), 
                  (-1, -1), (-1, 1), (1, -1), (1, 1)]
    opponent = '1' if player == '0' else '0'
    n = len(grid)
    
    def is_valid(nx, ny):
        return 0 <= nx < n and 0 <= ny < n
    
    def flip_in_direction(dx, dy):
        nx, ny = x + dx, y + dy
        pieces_to_flip = []
        while is_valid(nx, ny) and grid[nx][ny] == opponent:
            pieces_to_flip.append((nx, ny))
            nx += dx
            ny += dy
        if is_valid(nx, ny) and grid[nx][ny] == player:
            for px, py in pieces_to_flip:
                grid[px][py] = player
    
    for dx, dy in directions:
        flip_in_direction(dx, dy)

# Initial grid
grid = [
    ['*', '*', '*', '*'],
    ['0', '0', '0', '*'],
    ['*', '0', '1', '*'],
    ['*', '*', '*', '*']
]

# Place the black piece at (2, 1)
grid[2][1] = '0'

# Flip the pieces according to the rules
flip_pieces(grid, 2, 1, '0')

# Flatten the grid to the required format
flattened_grid = [cell for row in grid for cell in row]
print(flattened_grid)