# Initial grid
grid = [
    [0, 4, 0, 0, 0],
    [0, 0, 0, 0, 0],
    [0, 32, 0, 0, 0],
    [0, 0, 0, 16, 8],
    [0, 0, 0, 16, 16]
]

# Initial position
position = (0, 4)

# Movement sequence
moves = 'lrruru'

# Function to move the number
def move_number(grid, position, moves):
    x, y = position
    for move in moves:
        if move == 'l' and y > 0:
            y -= 1
        elif move == 'r' and y < 4:
            y += 1
        elif move == 'u' and x > 0:
            x -= 1
        elif move == 'd' and x < 4:
            x += 1
    return grid

# Simulate the movement
result_grid = move_number(grid, position, moves)

# Print the resulting grid
print(result_grid)