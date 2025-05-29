# Initial position
position = 0

# Take 4 steps forward
position += 4

# Take 2 more steps forward
position += 2

# Turn around and take 6 steps backward
position -= 6

# Check if the final position is the starting point
print("Yes" if position == 0 else "No")