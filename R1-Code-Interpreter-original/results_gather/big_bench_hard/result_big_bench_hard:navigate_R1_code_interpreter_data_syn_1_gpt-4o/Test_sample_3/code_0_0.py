# Initial position
x, y = 0, 0

# Step 1: Take 1 step right (east)
x += 1

# Step 2: Take 3 steps left (north)
y += 3

# Step 3: Take 2 steps right (east)
x += 2

# Check if we return to the starting point
return_to_start = (x == 0 and y == 0)

print("Yes" if return_to_start else "No")