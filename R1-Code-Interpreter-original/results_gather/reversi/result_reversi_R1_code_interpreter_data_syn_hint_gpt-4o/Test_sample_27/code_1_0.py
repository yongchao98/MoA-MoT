# Initial grid setup
grid = [
    [1, 0, '*', '*'],
    ['*', 1, 0, '*'],
    ['*', 0, 1, '*'],
    ['*', '*', '*', '*']
]

# Function to flip pieces
def flip_pieces(grid, x, y, player):
    opponent = 1 if player == 0 else 0
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (1, 1), (-1, 1), (1, -1)]
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

# Round 1: Your move (black) at (1,2)
grid[0][1] = 0
flip_pieces(grid, 0, 1, 0)

# Round 2: Opponent's move (white) at (1,1)
grid[0][0] = 1
flip_pieces(grid, 0, 0, 1)

# Convert the grid to a single string
result = ','.join(str(cell) for row in grid for cell in row)
print(result)