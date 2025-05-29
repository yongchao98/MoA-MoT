# Initial position
x, y = 0, 0

# Movements
x -= 10  # 10 steps left (west)
x -= 4   # 4 steps forward (west)
y -= 1   # 1 step left (south)
y += 4   # 4 steps backward (north)

# Check if we return to the starting point
return_to_start = (x == 0 and y == 0)

print("Yes" if return_to_start else "No")