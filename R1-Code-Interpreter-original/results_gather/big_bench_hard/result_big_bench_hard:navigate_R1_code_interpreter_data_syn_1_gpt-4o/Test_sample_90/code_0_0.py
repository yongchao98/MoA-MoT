# Initial position
x, y = 0, 0

# Movements
y += 10  # 10 steps forward
x += 8   # 8 steps right
x += 8   # 8 steps right
x -= 9   # 9 steps left
x += 6   # 6 steps right
y += 2   # 2 steps forward
y -= 9   # 9 steps backward

# Check if final position is the starting point
final_position = (x, y)
starting_point = (0, 0)

# Output whether the final position is the starting point
print(final_position == starting_point)