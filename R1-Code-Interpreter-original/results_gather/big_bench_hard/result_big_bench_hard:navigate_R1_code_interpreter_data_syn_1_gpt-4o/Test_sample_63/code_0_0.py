# Initial position
x, y = 0, 0

# Movements
x += 1  # 1 step right
x -= 9  # 9 steps left
x += 9  # 9 steps right
y -= 1  # 1 step backward
x += 3  # 3 steps right
y -= 6  # 6 steps backward

# Check if we are back at the starting point
result = (x == 0 and y == 0)
print("Yes" if result else "No")