# Initial position
x, y = 0, 0

# Initial direction (facing north)
direction = 'N'

# Steps in the initial direction (north)
y += 1  # Take 1 step
y += 5  # Take 5 steps
y += 1  # Take 1 step
y += 6  # Take 6 steps

# Turn left (now facing west)
direction = 'W'

# Steps in the new direction (west)
x -= 7  # Take 7 steps
x -= 2  # Take 2 steps
x -= 4  # Take 4 steps

# Check if we return to the starting point
result = (x == 0 and y == 0)

print("Yes" if result else "No")