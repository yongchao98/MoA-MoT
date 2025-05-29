# Initial position
position = 0

# Take 4 steps forward
position += 4

# Take 2 steps backward
position -= 2

# Take 10 steps backward
position -= 10

# Check if the final position is the same as the starting position
result = (position == 0)

print("Yes" if result else "No")