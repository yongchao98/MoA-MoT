# Initial position
x, y = 0, 0

# Movements
x += 9  # 9 steps right
x += 6  # 6 steps right
y -= 10 # 10 steps backward
x -= 9  # 9 steps left
x -= 4  # 4 steps left

# Final position
final_position = (x, y)

# Check if final position is the same as the starting point
starting_point = (0, 0)
result = final_position == starting_point

print("Yes" if result else "No")