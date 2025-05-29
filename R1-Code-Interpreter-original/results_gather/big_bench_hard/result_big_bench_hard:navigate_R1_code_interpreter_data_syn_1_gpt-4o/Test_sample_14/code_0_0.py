# Initial position
position = 0

# Take 3 steps left
position -= 3

# Take 10 steps right
position += 10

# Take 7 steps left
position -= 7

# Check if we return to the starting point
result = "Yes" if position == 0 else "No"
print(result)