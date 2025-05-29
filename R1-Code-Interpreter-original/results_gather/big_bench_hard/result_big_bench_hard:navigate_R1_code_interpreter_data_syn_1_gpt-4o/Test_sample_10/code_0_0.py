# Initial position
position = 0

# Steps taken
steps = [8, 5, 5, 10, 5]

# Calculate final position
for step in steps:
    position += step

# Check if the final position is the starting point
is_starting_point = (position == 0)

print(is_starting_point)