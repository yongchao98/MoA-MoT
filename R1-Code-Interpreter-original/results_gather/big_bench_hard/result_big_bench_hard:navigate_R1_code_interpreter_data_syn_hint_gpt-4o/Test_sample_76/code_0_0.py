# Initial position
x, y = 0, 0

# Movements
# Assume initial direction is north (positive y direction)
y += 1  # Take 1 step
y += 5  # Take 5 steps
y += 1  # Take 1 step
y += 6  # Take 6 steps

# Turn left, now facing west (negative x direction)
x -= 7  # Take 7 steps
x -= 2  # Take 2 steps
x -= 4  # Take 4 steps

# Check if we are back at the starting point
is_at_starting_point = (x == 0) and (y == 0)

print(is_at_starting_point)