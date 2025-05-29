# Initial position
x, y = 0, 0

# Step 1: Take 3 steps up
y += 3

# Step 2: Turn around and take 3 steps down
y -= 3

# Step 3: Turn right and take 3 steps
x += 3

# Check if the final position is the same as the starting point
final_position = (x, y)
starting_point = (0, 0)

# Output the result
print("Yes" if final_position == starting_point else "No")