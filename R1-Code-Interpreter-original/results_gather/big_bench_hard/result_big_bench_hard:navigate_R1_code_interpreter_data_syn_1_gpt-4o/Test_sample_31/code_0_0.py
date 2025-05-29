# Initial position
x, y = 0, 0

# Take 4 steps forward
y += 4

# Take 5 steps right
x += 5

# Take 7 steps left
x -= 7

# Check if returned to starting point
returned_to_start = (x == 0 and y == 0)

# Output the result
print("Yes" if returned_to_start else "No")