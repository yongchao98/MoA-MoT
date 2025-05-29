# Initial position
x, y = 0, 0

# Movements
movements = [
    (0, 7),  # 7 steps north
    (0, 7),  # 7 steps north
    (0, 10), # 10 steps north
    (0, 1),  # 1 step north
    (-7, 0)  # 7 steps west after turning left
]

# Calculate final position
for dx, dy in movements:
    x += dx
    y += dy

# Check if final position is the starting point
is_starting_point = (x == 0 and y == 0)
print(is_starting_point)