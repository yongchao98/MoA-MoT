# Initial position
x, y = 0, 0

# Movements
movements = [
    (0, -1),  # 1 step backward
    (5, 0),   # 5 steps right
    (-3, 0),  # 3 steps left
    (0, -1),  # 1 step backward
    (0, -5),  # 5 steps backward
    (0, 4),   # 4 steps forward
    (0, 3),   # 3 steps forward
    (-2, 0)   # 2 steps left
]

# Apply movements
for dx, dy in movements:
    x += dx
    y += dy

# Check if we return to the starting point
result = (x == 0 and y == 0)

print("Yes" if result else "No")