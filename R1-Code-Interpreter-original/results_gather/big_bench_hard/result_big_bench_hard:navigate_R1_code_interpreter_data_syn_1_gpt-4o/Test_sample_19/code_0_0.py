# Initial position
x, y = 0, 0

# Step 1: Take 8 steps forward (north)
y += 8

# Step 2: Turn left (now facing west) and take 9 steps
x -= 9

# Step 3: Take 7 steps forward (still facing west)
x -= 7

# Final position
print((x, y))