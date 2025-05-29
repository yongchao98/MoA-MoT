# Initial position and direction
x, y = 0, 0
direction = 'N'

# Define movements for each direction
movements = {
    'N': (0, 1),
    'E': (1, 0),
    'S': (0, -1),
    'W': (-1, 0)
}

# Define left turn changes
left_turn = {
    'N': 'W',
    'W': 'S',
    'S': 'E',
    'E': 'N'
}

# Steps and turns
steps = [3, 6, 6, 8, 1, 3, 7, 8]
turns = [4]  # Turn left after the 4th step

# Execute steps
for i, step in enumerate(steps):
    dx, dy = movements[direction]
    x += dx * step
    y += dy * step
    if i in turns:
        direction = left_turn[direction]

# Check if returned to starting point
returned_to_start = (x == 0 and y == 0)
print(returned_to_start)