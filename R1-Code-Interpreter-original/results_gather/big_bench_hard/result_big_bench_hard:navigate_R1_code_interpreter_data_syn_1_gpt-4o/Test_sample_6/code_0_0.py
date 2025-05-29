# Starting position
x, y = 0, 0

# Instructions
instructions = [
    ("right", 6),
    ("forward", 1),
    ("left", 10),
    ("forward", 8),
    ("backward", 9),
    ("right", 4)
]

# Process each instruction
for direction, steps in instructions:
    if direction == "right":
        x += steps
    elif direction == "left":
        x -= steps
    elif direction == "forward":
        y += steps
    elif direction == "backward":
        y -= steps

# Check if the final position is the starting point
is_starting_point = (x == 0 and y == 0)

print("Yes" if is_starting_point else "No")