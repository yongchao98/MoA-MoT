# Initial position
x, y = 0, 0

# Take 2 steps right
x += 2

# Take 2 steps right
x += 2

# Take 4 steps left
x -= 4

# Check if we are back to the starting point
result = (x == 0 and y == 0)
print("Yes" if result else "No")