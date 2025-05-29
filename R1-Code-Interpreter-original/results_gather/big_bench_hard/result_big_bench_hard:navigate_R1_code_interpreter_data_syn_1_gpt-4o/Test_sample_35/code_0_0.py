# Initial position
position = 0

# List of movements: positive for forward, negative for backward
movements = [7, -4, -4, 5, 7, -10, -1]

# Calculate the final position
for move in movements:
    position += move

# Output the final position
print(position)