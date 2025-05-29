# Initial position
x, y = 0, 0

# Take 6 steps in the positive x-direction
x += 6

# Turn right twice (facing negative x-direction) and take 7 steps
x -= 7

# Turn around (facing positive x-direction) and take 1 step
x += 1

# Check if we are back at the starting point
result = (x == 0 and y == 0)
print(result)