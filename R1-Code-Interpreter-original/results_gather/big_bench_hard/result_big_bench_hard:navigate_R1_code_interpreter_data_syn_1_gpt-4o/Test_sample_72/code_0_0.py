# Initialize position
position = 0

# Simulate the steps
position += 5  # Take 5 steps
position += 5  # Take 5 steps
position += 2  # Take 2 steps
position -= 8  # Turn around and take 8 steps
position -= 4  # Take 4 steps

# Check if returned to the starting point
print(position == 0)