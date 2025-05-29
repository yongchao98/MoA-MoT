# Initial position
x, y = 0, 0

# Move 9 steps left (negative y-direction)
y -= 9

# Move 4 steps backward (negative x-direction)
x -= 4

# Check if the final position is the starting point
is_starting_point = (x == 0 and y == 0)

print("Yes" if is_starting_point else "No")