# Initial position and direction
position = [0, 0]
direction = 'north'

# Function to move in the current direction
def move(steps):
    global position, direction
    if direction == 'north':
        position[1] += steps
    elif direction == 'east':
        position[0] += steps
    elif direction == 'south':
        position[1] -= steps
    elif direction == 'west':
        position[0] -= steps

# Function to turn right
def turn_right():
    global direction
    if direction == 'north':
        direction = 'east'
    elif direction == 'east':
        direction = 'south'
    elif direction == 'south':
        direction = 'west'
    elif direction == 'west':
        direction = 'north'

# Function to turn around
def turn_around():
    global direction
    if direction == 'north':
        direction = 'south'
    elif direction == 'east':
        direction = 'west'
    elif direction == 'south':
        direction = 'north'
    elif direction == 'west':
        direction = 'east'

# Simulate the steps
move(1)
move(10)
move(4)
move(1)
move(10)
turn_right()
turn_around()

# Check if we are back at the starting point
is_at_start = (position == [0, 0])
print("Yes" if is_at_start else "No")