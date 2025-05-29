# Initial position
x, y = 0, 0

# Movements
y -= 10  # 10 steps backward
y -= 7   # 7 steps backward
x += 8   # 8 steps right
x += 6   # 6 steps right
x -= 3   # 3 steps left

# Check if we return to the starting point
result = (x == 0 and y == 0)
print("Yes" if result else "No")