# Initial position and direction
x, y = 0, 0
direction = 'N'

# Function to move in the current direction
def move(steps, direction, x, y):
    if direction == 'N':
        y += steps
    elif direction == 'S':
        y -= steps
    elif direction == 'E':
        x += steps
    elif direction == 'W':
        x -= steps
    return x, y

# Function to turn around
def turn_around(direction):
    if direction == 'N':
        return 'S'
    elif direction == 'S':
        return 'N'
    elif direction == 'E':
        return 'W'
    elif direction == 'W':
        return 'E'

# Function to turn right
def turn_right(direction):
    if direction == 'N':
        return 'E'
    elif direction == 'E':
        return 'S'
    elif direction == 'S':
        return 'W'
    elif direction == 'W':
        return 'N'

# Function to turn left
def turn_left(direction):
    if direction == 'N':
        return 'W'
    elif direction == 'W':
        return 'S'
    elif direction == 'S':
        return 'E'
    elif direction == 'E':
        return 'N'

# Execute the steps
x, y = move(2, direction, x, y)  # Take 2 steps
direction = turn_around(direction)  # Turn around
direction = turn_right(direction)  # Turn right
direction = turn_left(direction)  # Turn left
x, y = move(2, direction, x, y)  # Take 2 steps
direction = turn_around(direction)  # Turn around
direction = turn_around(direction)  # Turn around

# Check if back to starting point
result = (x == 0 and y == 0)
print(result)