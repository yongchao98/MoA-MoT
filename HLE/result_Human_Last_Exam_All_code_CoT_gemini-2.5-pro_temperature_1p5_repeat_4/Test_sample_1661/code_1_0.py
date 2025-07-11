import sys

def solve():
    """
    Determines if a body is in equilibrium based on a series of force vectors.
    """
    try:
        # Read the number of force vectors.
        n = int(sys.stdin.readline())
    except (ValueError, IndexError):
        # If n is not a valid number, we can assume it's not in equilibrium.
        print("NO")
        return

    # Initialize the sum of force components for each axis.
    sum_x = 0
    sum_y = 0
    sum_z = 0

    # Loop through each force vector.
    for _ in range(n):
        try:
            # Read the three components of the force vector.
            line = sys.stdin.readline().split()
            # If the line is empty, stop processing.
            if not line:
                break
            
            x, y, z = map(int, line)

            # Add the components to the total sums.
            sum_x += x
            sum_y += y
            sum_z += z
        except (ValueError, IndexError):
            # If a line is malformed, the body is not in a predictable state.
            # We set a sum to a non-zero value to ensure the output is "NO".
            sum_x = 1 
            break

    # A body is in equilibrium if the net force vector is (0, 0, 0).
    # This means the sum of components on each axis must be zero.
    if sum_x == 0 and sum_y == 0 and sum_z == 0:
        print("YES")
    else:
        print("NO")

solve()