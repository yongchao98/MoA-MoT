# Initial grid
grid = [
    [0, 16, 0, 32, 0],
    [0, 0, 0, 0, 0],
    [32, 2, 8, 32, 0],
    [0, 32, 0, 0, 0],
    [0, 0, 0, 0, 0]
]

# Initial position of the number to move
position = (4, 5)

# Movement sequence
sequence = 'rdrrl'

# Function to move the number in the grid
def move_number(grid, position, sequence):
    x, y = position
    for move in sequence:
        if move == 'r':
            y = min(y + 1, 4)
        elif move == 'd':
            x = min(x + 1, 4)
        elif move == 'l':
            y = max(y - 1, 0)
        elif move == 'u':
            x = max(x - 1, 0)
    return grid

# Apply the movement
final_grid = move_number(grid, position, sequence)

# Print the final grid
print(final_grid)