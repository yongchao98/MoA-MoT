# Initial position
x, y = 0, 0

# Take 7 steps north
y += 7

# Turn right (facing east), then right again (facing south)
# Take 6 steps south
y -= 6

# Final position
final_position = (x, y)

# Check if final position is the starting point
is_starting_point = final_position == (0, 0)

print(is_starting_point)