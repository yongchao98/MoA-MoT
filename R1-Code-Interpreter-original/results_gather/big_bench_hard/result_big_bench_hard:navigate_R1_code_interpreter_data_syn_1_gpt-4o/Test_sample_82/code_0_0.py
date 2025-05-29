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

# Execute the steps
move(8)  # Take 8 steps
move(1)  # Take 1 step
turn_left()  # Turn left
turn_left()  # Turn left
move(9)  # Take 9 steps

# Check if we are back at the starting point
result = (x == 0 and y == 0)
print("Yes" if result else "No")