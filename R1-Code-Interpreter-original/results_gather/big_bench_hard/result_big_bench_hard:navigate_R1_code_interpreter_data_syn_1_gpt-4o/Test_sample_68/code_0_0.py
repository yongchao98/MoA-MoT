# Initial position and direction
position = [0, 0]
direction = 'N'

# Function to move in the current direction
def move(steps):
    global position, direction
    if direction == 'N':
        position[1] += steps
    elif direction == 'S':
        position[1] -= steps
    elif direction == 'E':
        position[0] += steps
    elif direction == 'W':
        position[0] -= steps

# Function to turn left
def turn_left():
    global direction
    if direction == 'N':
        direction = 'W'
    elif direction == 'W':
        direction = 'S'
    elif direction == 'S':
        direction = 'E'
    elif direction == 'E':
        direction = 'N'

# Function to turn around
def turn_around():
    global direction
    if direction == 'N':
        direction = 'S'
    elif direction == 'S':
        direction = 'N'
    elif direction == 'E':
        direction = 'W'
    elif direction == 'W':
        direction = 'E'

# Execute the instructions
move(3)
turn_left()
turn_left()
turn_around()
move(5)
turn_around()
move(5)

# Check if we are back at the starting point
result = position == [0, 0]
print("Yes" if result else "No")