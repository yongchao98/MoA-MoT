# Initial position
x, y = 0, 0

# Apply the movements
y -= 5  # 5 steps backward
x -= 2  # 2 steps left
x += 5  # 5 steps right
x -= 3  # 3 steps left
y -= 5  # 5 steps backward
y += 10 # 10 steps forward

# Check if the final position is the starting point
is_starting_point = (x == 0 and y == 0)

print("Yes" if is_starting_point else "No")