# Initial position
position = 0

# Take 6 steps forward
position += 6

# Take 2 steps backward
position -= 2

# Take 4 steps backward
position -= 4

# Check if we return to the starting point
result = "Yes" if position == 0 else "No"
print(result)