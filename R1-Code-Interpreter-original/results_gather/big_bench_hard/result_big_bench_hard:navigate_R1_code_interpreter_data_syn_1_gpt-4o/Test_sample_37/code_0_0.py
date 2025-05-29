# Initial position and direction
position = [0, 0]
direction = 'N'

# Function to update position based on direction and steps
def move(position, direction, steps):
    if direction == 'N':
        position[1] += steps
    elif direction == 'E':
        position[0] += steps
    elif direction == 'S':
        position[1] -= steps
    elif direction == 'W':
        position[0] -= steps
    return position

# Function to update direction based on current direction and turn
def turn(direction, turn_direction):
    directions = ['N', 'E', 'S', 'W']
    idx = directions.index(direction)
    if turn_direction == 'right':
        direction = directions[(idx + 1) % 4]
    elif turn_direction == 'left':
        direction = directions[(idx - 1) % 4]
    return direction

# Instructions
instructions = [
    ('right', 6),
    ('right', 5),
    ('left', 7),
    ('left', 3),
    ('left', 7),
    ('right', 6)
]

# Process each instruction
for turn_direction, steps in instructions:
    direction = turn(direction, turn_direction)
    position = move(position, direction, steps)

# Check if the final position is the starting point
is_starting_point = (position == [0, 0])

print("Yes" if is_starting_point else "No")