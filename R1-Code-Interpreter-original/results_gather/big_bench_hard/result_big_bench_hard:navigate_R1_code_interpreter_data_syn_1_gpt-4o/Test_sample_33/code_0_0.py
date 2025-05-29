# Initial position and direction
x, y = 0, 0
direction = 'N'  # N = North, E = East, S = South, W = West

# Function to move in the current direction
def move(steps):
    global x, y, direction
    if direction == 'N':
        y += steps
    elif direction == 'E':
        x += steps
    elif direction == 'S':
        y -= steps
    elif direction == 'W':
        x -= steps

# Function to turn right
def turn_right():
    global direction
    if direction == 'N':
        direction = 'E'
    elif direction == 'E':
        direction = 'S'
    elif direction == 'S':
        direction = 'W'
    elif direction == 'W':
        direction = 'N'

# Function to turn around
def turn_around():
    global direction
    if direction == 'N':
        direction = 'S'
    elif direction == 'E':
        direction = 'W'
    elif direction == 'S':
        direction = 'N'
    elif direction == 'W':
        direction = 'E'

# Execute the steps
move(3)
move(5)
move(8)
turn_right()
turn_right()
move(5)
turn_around()
move(2)
move(4)

# Check if back to starting point
is_at_start = (x == 0 and y == 0)
print("Yes" if is_at_start else "No")