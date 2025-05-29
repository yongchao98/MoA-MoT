# Initial position and direction
position = [0, 0]
direction = 'forward'  # Assume forward is positive y-direction

# Define direction changes
left_turns = {
    'forward': 'left',
    'left': 'backward',
    'backward': 'right',
    'right': 'forward'
}

# Define movement vectors
movement_vectors = {
    'forward': (0, 1),
    'left': (-1, 0),
    'backward': (0, -1),
    'right': (1, 0)
}

# Steps to follow
steps = [
    ('left', 6),
    ('forward', 7),
    ('left', 8),
    ('left', 7),
    ('forward', 6),
    ('forward', 1),
    ('forward', 4)
]

# Simulate the movement
for turn, step_count in steps:
    if turn == 'left':
        direction = left_turns[direction]
    dx, dy = movement_vectors[direction]
    position[0] += dx * step_count
    position[1] += dy * step_count

# Check if we return to the starting point
result = "Yes" if position == [0, 0] else "No"
print(result)