# Initial position
x, y = 0, 0

# Movements
movements = [
    (0, 9),   # 9 steps forward
    (7, 0),   # 7 steps left
    (0, -5),  # 5 steps backward
    (-4, 0),  # 4 steps right
    (-3, 0),  # 3 steps right
    (0, -4)   # 4 steps backward
]

# Apply movements
for dx, dy in movements:
    x += dx
    y += dy

# Check if we are back at the starting point
result = "Yes" if (x, y) == (0, 0) else "No"
print(result)