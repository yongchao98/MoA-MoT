# Initial position
position = 0

# Take 1 step
position += 1

# Take 5 steps
position += 5

# Turn around (change direction)
direction = -1

# Turn around (change direction again)
direction *= -1

# Turn around (change direction once more)
direction *= -1

# Take 6 steps in the current direction (backward)
position += 6 * direction

# Check if we are back at the starting point
result = "Yes" if position == 0 else "No"
print(result)