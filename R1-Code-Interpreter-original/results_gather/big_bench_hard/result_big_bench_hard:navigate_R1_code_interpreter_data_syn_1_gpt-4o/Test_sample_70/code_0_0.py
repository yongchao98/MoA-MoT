# Initial position
position = 0

# Take 8 steps backward
position -= 8

# Take 2 steps backward
position -= 2

# Take 10 steps forward
position += 10

# Check if the final position is the starting point
result = "Yes" if position == 0 else "No"
print(result)