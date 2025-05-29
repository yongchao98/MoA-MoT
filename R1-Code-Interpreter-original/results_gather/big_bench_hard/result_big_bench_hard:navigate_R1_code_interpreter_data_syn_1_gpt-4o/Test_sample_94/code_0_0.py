# Initial position and direction
position = [0, 0]
direction = 'N'

# Direction vectors
directions = {
    'N': (0, 1),
    'W': (-1, 0),
    'S': (0, -1),
    'E': (1, 0)
}

# Turn left function
def turn_left(current_direction):
    if current_direction == 'N':
        return 'W'
    elif current_direction == 'W':
        return 'S'
    elif current_direction == 'S':
        return 'E'
    elif current_direction == 'E':
        return 'N'

# Execute instructions
direction = turn_left(direction)  # Turn left
direction = turn_left(direction)  # Turn left
position[0] += directions[direction][0] * 10  # Take 10 steps
position[1] += directions[direction][1] * 10
position[0] += directions[direction][0] * 6  # Take 6 steps
position[1] += directions[direction][1] * 6
direction = turn_left(direction)  # Turn left
position[0] += directions[direction][0] * 3  # Take 3 steps
position[1] += directions[direction][1] * 3
position[0] += directions[direction][0] * 3  # Take 3 steps
position[1] += directions[direction][1] * 3

# Check if back to start
is_back_to_start = (position == [0, 0])
print("Yes" if is_back_to_start else "No")