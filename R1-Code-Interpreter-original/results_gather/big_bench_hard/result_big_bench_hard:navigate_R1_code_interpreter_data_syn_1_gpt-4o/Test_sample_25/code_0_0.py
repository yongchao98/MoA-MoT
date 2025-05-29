# Initial position and direction
x, y = 0, 0
direction = 'N'  # N = North, E = East, S = South, W = West

# Define movements based on direction
movements = {
    'N': (0, 1),
    'E': (1, 0),
    'S': (0, -1),
    'W': (-1, 0)
}

# Define turns
left_turn = {'N': 'W', 'W': 'S', 'S': 'E', 'E': 'N'}
right_turn = {'N': 'E', 'E': 'S', 'S': 'W', 'W': 'N'}

# Execute the steps
steps = [
    (9, None),  # Take 9 steps
    (6, None),  # Take 6 steps
    (9, None),  # Take 9 steps
    (6, None),  # Take 6 steps
    (2, None),  # Take 2 steps
    (None, 'L'),  # Turn left
    (6, None),  # Take 6 steps
    (None, 'A'),  # Turn around
    (8, None)   # Take 8 steps
]

for step, turn in steps:
    if step is not None:
        dx, dy = movements[direction]
        x += dx * step
        y += dy * step
    if turn == 'L':
        direction = left_turn[direction]
    elif turn == 'A':
        direction = right_turn[right_turn[direction]]

# Check if returned to starting point
returned_to_start = (x == 0 and y == 0)
print("Returned to starting point:", returned_to_start)