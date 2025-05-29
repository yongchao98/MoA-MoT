# Initial position
position = 0

# Take 7 steps left (negative direction)
position -= 7

# Take 2 steps backward (positive direction)
position += 2

# Take 7 steps backward (positive direction)
position += 7

# Take 7 steps backward (positive direction)
position += 7

# Take 3 steps forward (negative direction)
position -= 3

# Check if the final position is the starting point
is_starting_point = (position == 0)

print("Yes" if is_starting_point else "No")