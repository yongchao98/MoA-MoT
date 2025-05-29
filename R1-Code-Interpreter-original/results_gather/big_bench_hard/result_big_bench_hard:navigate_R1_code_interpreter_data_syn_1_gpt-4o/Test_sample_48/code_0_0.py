# Initial position and direction
position = [0, 0]  # [x, y] coordinates
direction = 'forward'  # Initial direction

# Define movements
movements = [
    ('right', 2),
    ('right', 4),
    ('left', 6)
]

# Define direction changes
direction_changes = {
    'forward': {'right': 'right', 'left': 'left'},
    'right': {'right': 'backward', 'left': 'forward'},
    'backward': {'right': 'left', 'left': 'right'},
    'left': {'right': 'forward', 'left': 'backward'}
}

# Define movement vectors
movement_vectors = {
    'forward': (0, 1),
    'right': (1, 0),
    'backward': (0, -1),
    'left': (-1, 0)
}

# Simulate movements
for turn, steps in movements:
    direction = direction_changes[direction][turn]
    move_vector = movement_vectors[direction]
    position[0] += move_vector[0] * steps
    position[1] += move_vector[1] * steps

# Check if returned to starting point
returned_to_start = position == [0, 0]
print("Yes" if returned_to_start else "No")