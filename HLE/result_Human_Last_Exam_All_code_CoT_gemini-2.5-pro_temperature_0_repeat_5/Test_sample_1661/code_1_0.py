import sys

def solve():
    """
    Determines if a body is in equilibrium based on the sum of force vectors.
    """
    try:
        # Read the number of force vectors.
        n = int(sys.stdin.readline())
    except (IOError, ValueError):
        # Handle cases with no input or invalid input.
        print("NO")
        return

    # Initialize sums for each coordinate.
    sum_x = 0
    sum_y = 0
    sum_z = 0

    # Loop n times to read each vector and add it to the sum.
    for _ in range(n):
        try:
            # Read the three components of the force vector.
            line = sys.stdin.readline()
            # Split the line and convert components to integers.
            components = list(map(int, line.split()))
            if len(components) == 3:
                x, y, z = components
                # Add components to the total sums.
                sum_x += x
                sum_y += y
                sum_z += z
            else:
                # Handle malformed lines
                print("NO")
                return
        except (IOError, ValueError):
            # Handle cases with incomplete or invalid input lines.
            print("NO")
            return

    # A body is in equilibrium if the sum of all force vectors is the zero vector.
    # This means the sum of components along each axis must be zero.
    if sum_x == 0 and sum_y == 0 and sum_z == 0:
        print("YES")
    else:
        print("NO")

solve()