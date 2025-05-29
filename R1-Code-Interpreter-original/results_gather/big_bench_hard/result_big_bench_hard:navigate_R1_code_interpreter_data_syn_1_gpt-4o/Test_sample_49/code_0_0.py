# Initial position
x, y = 0, 0

# Movements
x += 7  # 7 steps right (east)
x += 1  # 1 step forward (east)
y += 3  # 3 steps left (north)
x -= 2  # 2 steps left (west)

# Check if we return to the starting point
return_to_start = (x == 0 and y == 0)
print("Yes" if return_to_start else "No")