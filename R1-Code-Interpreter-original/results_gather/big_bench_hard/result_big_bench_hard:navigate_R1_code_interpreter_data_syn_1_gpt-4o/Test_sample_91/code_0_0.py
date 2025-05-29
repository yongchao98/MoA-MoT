# Initial position
x, y = 0, 0

# Turn left (from north to west)
# Move 9 steps west
x -= 9

# Move another 9 steps west
x -= 9

# Check if we are back at the starting point
result = (x == 0 and y == 0)
print("Yes" if result else "No")