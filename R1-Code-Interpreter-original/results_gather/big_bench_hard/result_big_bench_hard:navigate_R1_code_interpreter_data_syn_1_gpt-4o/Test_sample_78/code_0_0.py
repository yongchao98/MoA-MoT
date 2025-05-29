# Initial position
x, y = 0, 0

# Take 5 steps right
x += 5

# Take 8 steps backward
y -= 8

# Take 1 step right
x += 1

# Take 7 steps backward
y -= 7

# Check if we return to the starting point
result = (x == 0 and y == 0)
print("Yes" if result else "No")