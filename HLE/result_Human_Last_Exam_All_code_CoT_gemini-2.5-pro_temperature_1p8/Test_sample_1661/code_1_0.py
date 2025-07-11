import sys

def check_equilibrium():
    """
    Reads force vectors from standard input and determines if the net force is zero.
    """
    try:
        # Read the number of force vectors.
        n = int(sys.stdin.readline())
    except (ValueError, IndexError):
        # Handle cases with no or invalid input for n.
        print("NO")
        return

    # Initialize sums for each coordinate.
    sum_x, sum_y, sum_z = 0, 0, 0

    # Loop n times to read each force vector.
    for _ in range(n):
        try:
            # Read and parse the coordinates of a force vector.
            xi, yi, zi = map(int, sys.stdin.readline().split())
            
            # Add the components to the total sums.
            sum_x += xi
            sum_y += yi
            sum_z += zi
        except (ValueError, IndexError):
            # Handle malformed lines.
            print("NO")
            return
            
    # A body is in equilibrium if the sum of forces on each axis is zero.
    if sum_x == 0 and sum_y == 0 and sum_z == 0:
        print("YES")
    else:
        print("NO")

check_equilibrium()