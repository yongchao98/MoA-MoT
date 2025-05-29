# Initial position
x, y = 0, 0

# Step 1: Take 3 steps right (facing north, right is east)
x += 3

# Step 2: Take 5 steps right (facing east, right is south)
y -= 5

# Step 3: Take 2 steps backward (facing south, backward is north)
y += 2

# Final position
final_position = (x, y)

# Check if the final position is the starting point
is_starting_point = final_position == (0, 0)

print("Final Position:", final_position)
print("Return to Starting Point:", is_starting_point)