import sys

# Step 1: Read the number of force vectors from input.
try:
    n = int(sys.stdin.readline())
except (ValueError, IndexError):
    # If the first line is not a valid number, it's an invalid input.
    # The net force is undefined, so it's not in a state of equilibrium.
    print("NO")
    sys.exit()

# Step 2: Initialize sums of force components to zero.
sum_x = 0
sum_y = 0
sum_z = 0

# Step 3 & 4: Loop n times to read each vector and add its components to the sums.
for _ in range(n):
    try:
        # Read the three integer components from a single line.
        x, y, z = map(int, sys.stdin.readline().split())
        
        # Add the components to their respective sums.
        sum_x += x
        sum_y += y
        sum_z += z
    except (ValueError, IndexError):
        # Assuming valid input format as per the problem description.
        # If an invalid line is encountered, we can treat the sums as incorrect.
        # For simplicity, we can pass, but a robust solution might exit.
        pass

# Step 5 & 6: Check if the body is in equilibrium and print the result.
# The condition for equilibrium is that the net force is the zero vector,
# which means the sum of components along each axis must be zero.
if sum_x == 0 and sum_y == 0 and sum_z == 0:
    print("YES")
else:
    print("NO")