# Initial grid
grid = [
    [0, 16, 0, 32, 0],
    [0, 0, 0, 0, 0],
    [32, 2, 8, 32, 0],
    [0, 32, 0, 0, 0],
    [0, 0, 0, 0, 0]
]

# Function to move and combine numbers in the grid
def move_and_combine(grid, sequence):
    # Define the movement directions
    directions = {
        'r': (0, 1),
        'd': (1, 0),
        'l': (0, -1),
        'u': (-1, 0)
    }
    
    # Function to apply a single move
    def apply_move(grid, direction):
        dx, dy = directions[direction]
        new_grid = [[0] * 5 for _ in range(5)]
        for x in range(5):
            for y in range(5):
                if grid[x][y] != 0:
                    new_x, new_y = x, y
                    while 0 <= new_x + dx < 5 and 0 <= new_y + dy < 5 and new_grid[new_x + dx][new_y + dy] == 0:
                        new_x += dx
                        new_y += dy
                    if 0 <= new_x + dx < 5 and 0 <= new_y + dy < 5 and new_grid[new_x + dx][new_y + dy] == grid[x][y]:
                        new_grid[new_x + dx][new_y + dy] *= 2
                    else:
                        new_grid[new_x][new_y] = grid[x][y]
        return new_grid
    
    # Apply the sequence of moves
    for move in sequence:
        grid = apply_move(grid, move)
    
    return grid

# Apply the movement and combination
final_grid = move_and_combine(grid, 'rdrrl')

# Print the final grid
print(final_grid)