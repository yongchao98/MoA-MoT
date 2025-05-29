# Initial position and direction
position = [0, 0]
directions = [(0, 1), (-1, 0), (0, -1), (1, 0)]  # North, West, South, East
current_direction = 0  # Start facing North

# Instructions
instructions = [
    ("left", 7),
    ("left", 8),
    ("left", 1),
    ("left", 9),
    ("forward", 10),
    ("forward", 5)
]

# Process each instruction
for action, steps in instructions:
    if action == "left":
        current_direction = (current_direction + 1) % 4  # Turn left
    elif action == "forward":
        dx, dy = directions[current_direction]
        position[0] += dx * steps
        position[1] += dy * steps

# Check if the final position is the starting point
result = position == [0, 0]
print("Yes" if result else "No")