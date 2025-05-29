def print_grid(grid):
    for row in grid:
        print(''.join(row))

def apply_move(grid, move, player):
    n = len(grid)
    x, y = move
    grid[x][y] = player

    directions = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (-1, 1), (1, -1), (1, 1)]

    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        pieces_to_flip = []

        while 0 <= nx < n and 0 <= ny < n:
            if grid[nx][ny] == '*':
                break
            if grid[nx][ny] == player:
                for px, py in pieces_to_flip:
                    grid[px][py] = player
                break
            pieces_to_flip.append((nx, ny))
            nx += dx
            ny += dy

# Initial grid
grid = [
    ['*', '*', '*', '*'],
    ['0', '0', '0', '*'],
    ['1', '0', '0', '*'],
    ['*', '*', '0', '*']
]

# Moves: (row, column) are 0-indexed
moves = [(1, 0, '0'), (2, 0, '1'), (3, 2, '0')]

# Apply each move
for move in moves:
    x, y, player = move
    apply_move(grid, (x, y), player)

# Print the final grid
print_grid(grid)