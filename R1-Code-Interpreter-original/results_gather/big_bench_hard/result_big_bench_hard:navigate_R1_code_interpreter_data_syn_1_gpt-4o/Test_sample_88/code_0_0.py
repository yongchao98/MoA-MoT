# Initial position and direction
x, y = 0, 0
direction = 'N'  # N = North, E = East, S = South, W = West

# Function to turn left
def turn_left(current_direction):
    directions = ['N', 'W', 'S', 'E']
    return directions[(directions.index(current_direction) + 1) % 4]

# Function to turn right
def turn_right(current_direction):
    directions = ['N', 'E', 'S', 'W']
    return directions[(directions.index(current_direction) + 1) % 4]

# Function to move forward
def move_forward(x, y, direction, steps):
    if direction == 'N':
        y += steps
    elif direction == 'E':
        x += steps
    elif direction == 'S':
        y -= steps
    elif direction == 'W':
        x -= steps
    return x, y

# Follow the instructions
direction = turn_left(direction)
direction = turn_right(direction)
x, y = move_forward(x, y, direction, 8)
x, y = move_forward(x, y, direction, 4)
direction = turn_right(direction)
direction = turn_right(direction)
x, y = move_forward(x, y, direction, 10)
x, y = move_forward(x, y, direction, 1)
x, y = move_forward(x, y, direction, 1)

# Check if we are back at the starting point
result = (x == 0 and y == 0)
print("Yes" if result else "No")