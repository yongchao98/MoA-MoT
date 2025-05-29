# Initial grid setup
grid = [
    ['1', '0', '*', '0'],
    ['0', '0', '0', '*'],
    ['*', '0', '1', '*'],
    ['*', '*', '*', '*']
]

# Function to flip pieces
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

# Round 1: Black places at (2,1)
grid[1][0] = '0'
flip_pieces(grid, 1, 0, '0')

# Round 2: White places at (1,3)
grid[0][2] = '1'
flip_pieces(grid, 0, 2, '1')

# Round 3: Black places at (1,4)
grid[0][3] = '0'
flip_pieces(grid, 0, 3, '0')

# Flatten the grid to the required format
flattened_grid = [cell for row in grid for cell in row]
formatted_output = ','.join(flattened_grid)
print(formatted_output)