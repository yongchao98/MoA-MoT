# Initial position
x, y = 0, 0

# Take 10 steps left (west)
x -= 10

# Take 4 steps forward (west)
x -= 4

# Take 1 step left (south)
y -= 1

# Take 4 steps backward (north)
y += 4

# Check if the final position is the same as the starting point
result = (x == 0 and y == 0)
print("Yes" if result else "No")