def flip_pieces(grid, x, y, player):
    # Directions: right, left, down, up, down-right, down-left, up-right, up-left
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
    ['0', '0', '0', '*'],
    ['*', '*', '*', '*']
]

# Place the black piece at (3,4)
grid[2][3] = '0'

# Apply flipping rules
flip_pieces(grid, 2, 3, '0')

# Convert the grid to a single string representation
result = ','.join([cell for row in grid for cell in row])
print(result)