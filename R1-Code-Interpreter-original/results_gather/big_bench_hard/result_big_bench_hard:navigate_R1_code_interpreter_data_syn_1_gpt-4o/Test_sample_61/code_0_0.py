# Initial position
x, y = 0, 0

# Movements
movements = [
    (0, 10),  # 10 steps forward
    (10, 0),  # 10 steps right
    (0, -10), # 10 steps backward
    (0, 3),   # 3 steps forward
    (0, 7),   # 7 steps forward
    (-9, 0)   # 9 steps left
]

# Apply movements
for dx, dy in movements:
    x += dx
    y += dy

# Check if final position is the starting point
result = (x, y) == (0, 0)

print("Yes" if result else "No")