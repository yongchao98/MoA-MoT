# Initial position
position = 0

# Steps
steps = [6, 2, 6]

# Simulate the steps
for step in steps:
    position += step

# Check if we return to the starting point
print("Yes" if position == 0 else "No")