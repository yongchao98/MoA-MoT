# Initial position
x, y = 0, 0

# Move 8 steps forward
y += 8

# Move 9 steps right
x += 9

# Move 2 steps forward
y += 2

# Move 10 steps backward
y -= 10

# Move 9 steps left
x -= 9

# Check if we are back at the starting point
result = (x == 0 and y == 0)
print("Yes" if result else "No")