# Initial position
position = 0

# Take 7 steps forward
position += 7

# Take 3 steps backward
position -= 3

# Take 4 steps backward
position -= 4

# Check if we return to the starting point
result = "Yes" if position == 0 else "No"
print(result)