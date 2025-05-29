# Initial position and direction
position = [0, 0]
direction = 'north'

# Function to turn right
def turn_right(current_direction):
    directions = ['north', 'east', 'south', 'west']
    return directions[(directions.index(current_direction) + 1) % 4]

# Function to turn around
def turn_around(current_direction):
    directions = ['north', 'east', 'south', 'west']
    return directions[(directions.index(current_direction) + 2) % 4]

# Function to move steps in the current direction
def move_steps(position, direction, steps):
    if direction == 'north':
        position[1] += steps
    elif direction == 'east':
        position[0] += steps
    elif direction == 'south':
        position[1] -= steps
    elif direction == 'west':
        position[0] -= steps
    return position

# Execute the instructions
direction = turn_right(direction)  # Turn right
position = move_steps(position, direction, 9)  # Take 9 steps
direction = turn_around(direction)  # Turn around
position = move_steps(position, direction, 9)  # Take 9 steps
direction = turn_right(direction)  # Turn right

# Check if returned to the starting point
returned_to_start = (position == [0, 0])

print("Yes" if returned_to_start else "No")