# Initial position
x, y = 0, 0

# Movements
y += 8  # 8 steps forward
y -= 5  # 5 steps backward
x -= 5  # 5 steps left
y += 5  # 5 steps forward

# Check if the final position is the same as the starting point
final_position = (x, y)
starting_point = (0, 0)

# Output the result
print("Yes" if final_position == starting_point else "No")