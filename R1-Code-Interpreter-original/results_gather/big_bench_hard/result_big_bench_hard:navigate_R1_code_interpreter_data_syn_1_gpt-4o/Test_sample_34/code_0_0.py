# Initial position and direction
x, y = 0, 0
direction = 'N'

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

# Execute the instructions
move(7)
move(8)
turn_right()
turn_right()
move(1)
move(7)
move(6)
move(1)

# Check if we are back at the starting point
result = (x == 0 and y == 0)
print("Yes" if result else "No")