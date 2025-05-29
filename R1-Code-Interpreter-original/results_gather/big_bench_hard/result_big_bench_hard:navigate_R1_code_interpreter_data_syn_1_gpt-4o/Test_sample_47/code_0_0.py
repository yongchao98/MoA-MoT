# Initial position and direction
position = [0, 0]
direction = 'north'

# Define direction changes
direction_changes = {
    'north': {'left': 'west', 'right': 'east', 'around': 'south'},
    'south': {'left': 'east', 'right': 'west', 'around': 'north'},
    'east': {'left': 'north', 'right': 'south', 'around': 'west'},
    'west': {'left': 'south', 'right': 'north', 'around': 'east'}
}

# Define movement vectors
movement_vectors = {
    'north': (0, 1),
    'south': (0, -1),
    'east': (1, 0),
    'west': (-1, 0)
}

# Execute instructions
instructions = [
    ('turn', 'left'),
    ('turn', 'around'),
    ('step', 1),
    ('step', 10),
    ('turn', 'around'),
    ('step', 7),
    ('step', 4)
]

for action, value in instructions:
    if action == 'turn':
        direction = direction_changes[direction][value]
    elif action == 'step':
        move_x, move_y = movement_vectors[direction]
        position[0] += move_x * value
        position[1] += move_y * value

# Check if we return to the starting point
result = position == [0, 0]
print("Yes" if result else "No")