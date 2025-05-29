# Initial position
x, y = 0, 0

# Take 10 steps forward (initially along the positive x-axis)
x += 10

# Turn left (now facing positive y-axis)
# Turn left again (now facing negative x-axis)

# Take 10 steps in the negative x-axis direction
x -= 10

# Check if we are back at the starting point
result = (x, y) == (0, 0)
print("Yes" if result else "No")