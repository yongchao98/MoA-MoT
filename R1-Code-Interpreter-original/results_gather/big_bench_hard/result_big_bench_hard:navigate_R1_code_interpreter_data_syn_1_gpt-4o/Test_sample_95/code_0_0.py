# Initial position
position = 0

# Steps taken
steps = [2, 8, 8, 6]

# Calculate final position
for step in steps:
    position += step

# Check if the final position is the starting point
print("Yes" if position == 0 else "No")