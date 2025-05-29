# Initial position and direction
position = [0, 0]
direction = 'north'

# Define direction changes
right_turns = {'north': 'east', 'east': 'south', 'south': 'west', 'west': 'north'}
left_turns = {'north': 'west', 'west': 'south', 'south': 'east', 'east': 'north'}
opposite_directions = {'north': 'south', 'south': 'north', 'east': 'west', 'west': 'east'}

# Define movement vectors
movement_vectors = {'north': (0, 1), 'south': (0, -1), 'east': (1, 0), 'west': (-1, 0)}

# Turn right and take 3 steps
direction = right_turns[direction]
position[0] += 3 * movement_vectors[direction][0]
position[1] += 3 * movement_vectors[direction][1]

# Turn around and take 3 steps
direction = opposite_directions[direction]
position[0] += 3 * movement_vectors[direction][0]
position[1] += 3 * movement_vectors[direction][1]

# Turn left
direction = left_turns[direction]

# Check if we are back at the starting point
is_at_starting_point = (position == [0, 0])

print("Yes" if is_at_starting_point else "No")