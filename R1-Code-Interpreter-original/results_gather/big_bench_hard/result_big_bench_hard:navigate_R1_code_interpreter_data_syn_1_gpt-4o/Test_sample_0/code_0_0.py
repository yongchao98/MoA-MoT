# Initial position
x, y = 0, 0

# Movements
movements = [
    (0, -1),  # 1 step backward
    (-9, 0),  # 9 steps left
    (0, -2),  # 2 steps backward
    (0, 6),   # 6 steps forward
    (0, 4),   # 4 steps forward
    (0, -4),  # 4 steps backward
    (3, 0)    # 3 steps right
]

# Apply movements
for dx, dy in movements:
    x += dx
    y += dy

# Check if final position is the starting point
result = (x == 0 and y == 0)
print("Yes" if result else "No")