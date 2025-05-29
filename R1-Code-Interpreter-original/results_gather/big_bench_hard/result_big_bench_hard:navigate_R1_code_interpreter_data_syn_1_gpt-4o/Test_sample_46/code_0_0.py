# Initial position and direction
position = [0, 0]
direction = 'N'

# Define direction changes
right_turns = {'N': 'E', 'E': 'S', 'S': 'W', 'W': 'N'}
left_turns = {'N': 'W', 'W': 'S', 'S': 'E', 'E': 'N'}
move_steps = {'N': (0, 1), 'E': (1, 0), 'S': (0, -1), 'W': (-1, 0)}

# Instructions
instructions = [
    ('R', 0),  # Turn right
    ('S', 9),  # Take 9 steps
    ('S', 10), # Take 10 steps
    ('L', 0),  # Turn left
    ('S', 10), # Take 10 steps
    ('R', 0),  # Turn right
    ('R', 0),  # Turn around (right twice)
    ('S', 4)   # Take 4 steps
]

# Execute instructions
for action, value in instructions:
    if action == 'R':
        direction = right_turns[direction]
    elif action == 'L':
        direction = left_turns[direction]
    elif action == 'S':
        dx, dy = move_steps[direction]
        position[0] += dx * value
        position[1] += dy * value

# Check if returned to starting point
returned_to_start = (position == [0, 0])
print("Yes" if returned_to_start else "No")