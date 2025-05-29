# Initial position
x, y = 0, 0

# Movements
movements = [
    (0, -8),  # 8 steps backward
    (0, 9),   # 9 steps forward
    (0, 4),   # 4 steps forward
    (-2, 0),  # 2 steps left
    (0, 4),   # 4 steps forward
    (0, 6),   # 6 steps forward
    (0, -10), # 10 steps backward
    (8, 0),   # 8 steps right
    (-1, 0)   # 1 step left
]

# Apply movements
for dx, dy in movements:
    x += dx
    y += dy

# Check if we return to the starting point
result = (x == 0 and y == 0)
print("Yes" if result else "No")