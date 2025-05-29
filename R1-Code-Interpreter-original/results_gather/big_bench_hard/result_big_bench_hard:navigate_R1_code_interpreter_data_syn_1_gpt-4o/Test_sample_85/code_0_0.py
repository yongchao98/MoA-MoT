# Initial position
x, y = 0, 0

# Movements
y += 1  # 1 step forward
y += 10 # 10 steps forward
y -= 6  # 6 steps backward
y -= 7  # 7 steps backward
x -= 1  # 1 step left
x += 10 # 10 steps right

# Check if final position is the same as the starting point
final_position = (x, y)
starting_point = (0, 0)

# Output the result
print(final_position == starting_point)