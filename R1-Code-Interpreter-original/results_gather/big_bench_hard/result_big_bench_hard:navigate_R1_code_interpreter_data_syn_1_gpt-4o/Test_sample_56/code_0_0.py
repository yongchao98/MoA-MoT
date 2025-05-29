# Initial position and direction
position = [0, 0]
direction = 'N'

# Define movements for each direction
movements = {
    'N': (0, 1),
    'E': (1, 0),
    'S': (0, -1),
    'W': (-1, 0)
}

# Define left turn mapping
left_turn = {
    'N': 'W',
    'W': 'S',
    'S': 'E',
    'E': 'N'
}

# Execute the steps
# Take 9 steps
for _ in range(9):
    move = movements[direction]
    position[0] += move[0]
    position[1] += move[1]

# Turn left twice
direction = left_turn[direction]
direction = left_turn[direction]

# Take 5 steps
for _ in range(5):
    move = movements[direction]
    position[0] += move[0]
    position[1] += move[1]

# Take 4 steps
for _ in range(4):
    move = movements[direction]
    position[0] += move[0]
    position[1] += move[1]

# Turn left twice
direction = left_turn[direction]
direction = left_turn[direction]

# Check if we are back at the starting point
result = (position == [0, 0])
print("Yes" if result else "No")