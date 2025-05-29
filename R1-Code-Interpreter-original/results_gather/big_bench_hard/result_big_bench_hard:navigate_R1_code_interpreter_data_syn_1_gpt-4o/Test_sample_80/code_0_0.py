# Initial position and direction
position = [0, 0]
direction = 'N'

# Function to turn around
def turn_around(current_direction):
    if current_direction == 'N':
        return 'S'
    elif current_direction == 'S':
        return 'N'
    elif current_direction == 'E':
        return 'W'
    elif current_direction == 'W':
        return 'E'

# Function to turn left
def turn_left(current_direction):
    if current_direction == 'N':
        return 'W'
    elif current_direction == 'S':
        return 'E'
    elif current_direction == 'E':
        return 'N'
    elif current_direction == 'W':
        return 'S'

# Function to take steps
def take_steps(current_position, current_direction, steps):
    if current_direction == 'N':
        current_position[1] += steps
    elif current_direction == 'S':
        current_position[1] -= steps
    elif current_direction == 'E':
        current_position[0] += steps
    elif current_direction == 'W':
        current_position[0] -= steps
    return current_position

# Execute the instructions
direction = turn_around(direction)
direction = turn_left(direction)
position = take_steps(position, direction, 1)
position = take_steps(position, direction, 2)
direction = turn_around(direction)
direction = turn_around(direction)
direction = turn_around(direction)
position = take_steps(position, direction, 3)

# Check if we are back at the starting point
result = (position == [0, 0])
print(result)