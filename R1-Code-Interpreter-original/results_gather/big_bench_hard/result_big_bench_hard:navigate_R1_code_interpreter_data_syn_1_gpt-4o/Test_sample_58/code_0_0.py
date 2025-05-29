# Initial position and direction
x, y = 0, 0
direction = 'north'

# Define movements for each direction
movements = {
    'north': (0, 1),
    'east': (1, 0),
    'south': (0, -1),
    'west': (-1, 0)
}

# Function to move in the current direction
def move(steps):
    global x, y, direction
    dx, dy = movements[direction]
    x += dx * steps
    y += dy * steps

# Function to turn right
def turn_right():
    global direction
    directions = ['north', 'east', 'south', 'west']
    current_index = directions.index(direction)
    direction = directions[(current_index + 1) % 4]

# Execute the steps
move(8)  # Take 8 steps
move(1)  # Take 1 step
move(8)  # Take 8 steps
turn_right()  # Turn right

# Check if we are back at the starting point
result = (x == 0 and y == 0)
print("Yes" if result else "No")