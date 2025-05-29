# Initial position
x, y = 0, 0

# Movements
y -= 1  # 1 step backward
y += 5  # 5 steps forward
y -= 10 # 10 steps backward
y += 9  # 9 steps forward
x += 7  # 7 steps right

# Check if final position is the same as the starting point
final_position = (x, y)
starting_point = (0, 0)

# Output the result
print(final_position == starting_point)