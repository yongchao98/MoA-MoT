# Initialize position and direction
position = [0, 0]  # [x, y] where x is East-West, y is North-South
direction = 'N'  # Starting facing North

# Define direction changes
right_turns = {'N': 'E', 'E': 'S', 'S': 'W', 'W': 'N'}
left_turns = {'N': 'W', 'W': 'S', 'S': 'E', 'E': 'N'}

# Define movement vectors
movement = {'N': (0, 1), 'E': (1, 0), 'S': (0, -1), 'W': (-1, 0)}

# List of instructions
instructions = [
    (6, 'right'),
    (5, 'right'),
    (7, 'left'),
    (3, 'left'),
    (7, 'left'),
    (6, 'right')
]

# Process each instruction
for steps, turn in instructions:
    # Update direction
    if turn == 'right':
        direction = right_turns[direction]
    elif turn == 'left':
        direction = left_turns[direction]
    
    # Update position
    dx, dy = movement[direction]
    position[0] += dx * steps
    position[1] += dy * steps

# Check if returned to starting point
returned_to_start = position == [0, 0]
print("Yes" if returned_to_start else "No")