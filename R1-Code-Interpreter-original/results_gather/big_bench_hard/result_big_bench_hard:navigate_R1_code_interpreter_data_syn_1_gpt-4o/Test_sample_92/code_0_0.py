# Initial position and direction
x, y = 0, 0
direction = 'N'  # N for North, E for East, S for South, W for West

# Function to update position based on direction and steps
def move(x, y, direction, steps):
    if direction == 'N':
        y += steps
    elif direction == 'E':
        x += steps
    elif direction == 'S':
        y -= steps
    elif direction == 'W':
        x -= steps
    return x, y

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

# Execute the movements
# 1. Take 6 steps left
direction = turn_left(direction)
x, y = move(x, y, direction, 6)

# 2. Take 3 steps backward
# Backward is opposite to the current direction
if direction == 'N':
    x, y = move(x, y, 'S', 3)
elif direction == 'E':
    x, y = move(x, y, 'W', 3)
elif direction == 'S':
    x, y = move(x, y, 'N', 3)
elif direction == 'W':
    x, y = move(x, y, 'E', 3)

# 3. Take 3 steps left
direction = turn_left(direction)
x, y = move(x, y, direction, 3)

# Check if we are back at the starting point
result = (x == 0 and y == 0)
print("Yes" if result else "No")