# Initial matrix
matrix = [
    [0, 2, 0, 4, 0],
    [0, 16, 0, 0, 0],
    [0, 0, 0, 0, 0],
    [32, 0, 8, 0, 0],
    [0, 0, 16, 0, 0]
]

# Starting position
position = (1, 3)

# Movement sequence
moves = 'rurlrr'

# Define movement directions
directions = {
    'u': (-1, 0),
    'd': (1, 0),
    'l': (0, -1),
    'r': (0, 1)
}

# Simulate the movement
for move in moves:
    dx, dy = directions[move]
    new_x = position[0] + dx
    new_y = position[1] + dy
    # Check bounds
    if 0 <= new_x < 5 and 0 <= new_y < 5:
        position = (new_x, new_y)

# Print the final matrix
print(matrix)