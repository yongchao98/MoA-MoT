# Initial position
position = 0

# Take 10 steps forward
position += 10

# Take 4 steps backward
position -= 4

# Take 6 steps forward
position += 6

# Check if we return to the starting point
result = "Yes" if position == 0 else "No"
print(result)