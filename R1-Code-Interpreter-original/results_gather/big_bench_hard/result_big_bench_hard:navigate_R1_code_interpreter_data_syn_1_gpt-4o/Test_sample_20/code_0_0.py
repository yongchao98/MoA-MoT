# Initial position
x, y = 0, 0

# Take 6 steps backward
y -= 6

# Take 6 steps forward
y += 6

# Take 6 steps right
x += 6

# Check if the final position is the starting point
is_starting_point = (x == 0 and y == 0)

print("Yes" if is_starting_point else "No")