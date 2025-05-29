# Initial position and direction
position = [0, 0]
direction = 'north'

# Define direction changes
directions = ['north', 'west', 'south', 'east']
direction_index = directions.index(direction)

# Turn left
direction_index = (direction_index + 1) % 4
direction = directions[direction_index]

# Take 4 steps
if direction == 'north':
    position[1] += 4
elif direction == 'south':
    position[1] -= 4
elif direction == 'east':
    position[0] += 4
elif direction == 'west':
    position[0] -= 4

# Turn around
direction_index = (direction_index + 2) % 4
direction = directions[direction_index]

# Take 4 steps
if direction == 'north':
    position[1] += 4
elif direction == 'south':
    position[1] -= 4
elif direction == 'east':
    position[0] += 4
elif direction == 'west':
    position[0] -= 4

# Check if back to starting point
result = (position == [0, 0])
print("Yes" if result else "No")