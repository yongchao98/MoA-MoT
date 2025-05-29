# Initial position and direction
position = [0, 0]
direction = 'N'

# Define movements for each direction
movements = {
    'N': (0, 1),
    'S': (0, -1),
    'E': (1, 0),
    'W': (-1, 0)
}

# Define direction changes
right_turn = {'N': 'E', 'E': 'S', 'S': 'W', 'W': 'N'}
left_turn = {v: k for k, v in right_turn.items()}

# Execute the instructions
# Take 3 steps
for _ in range(3):
    position[0] += movements[direction][0]
    position[1] += movements[direction][1]

# Turn around
direction = left_turn[left_turn[direction]]

# Take 5 steps
for _ in range(5):
    position[0] += movements[direction][0]
    position[1] += movements[direction][1]

# Turn right
direction = right_turn[direction]

# Turn right
direction = right_turn[direction]

# Take 1 step
for _ in range(1):
    position[0] += movements[direction][0]
    position[1] += movements[direction][1]

# Take 1 step
for _ in range(1):
    position[0] += movements[direction][0]
    position[1] += movements[direction][1]

# Check if returned to starting point
returned_to_start = (position == [0, 0])
print("Yes" if returned_to_start else "No")