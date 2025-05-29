# Initial position
x, y = 0, 0

# Movements
movements = [
    (-9, 0),  # 9 steps left
    (0, 7),   # 7 steps forward
    (9, 0),   # 9 steps right
    (-7, 0),  # 7 steps left
    (0, -7),  # 7 steps backward
    (8, 0),   # 8 steps right
    (-10, 0), # 10 steps left
    (0, -1),  # 1 step backward
    (0, 1)    # 1 step forward
]

# Apply movements
for dx, dy in movements:
    x += dx
    y += dy

# Check if back to starting point
result = (x == 0 and y == 0)
print("Yes" if result else "No")