import sys

def solve():
    """
    This function determines if a body is in equilibrium by checking if the
    vector sum of all forces is the zero vector.
    """
    try:
        # Read the number of force vectors.
        n = int(sys.stdin.readline())
    except (ValueError, IndexError):
        # Handle cases with no or invalid input for n.
        print("NO")
        return

    # Initialize sums for each coordinate to zero.
    sum_x, sum_y, sum_z = 0, 0, 0

    # Loop n times to read all force vectors.
    for _ in range(n):
        try:
            # Read a line, split it, and convert components to integers.
            line = sys.stdin.readline()
            if not line:
                # Handle unexpected end of input.
                continue
            xi, yi, zi = map(int, line.split())

            # Add the components to the total sums.
            sum_x += xi
            sum_y += yi
            sum_z += zi
        except (ValueError, IndexError):
            # Handle malformed lines of input.
            continue

    # A body is in equilibrium if the net force along each axis is zero.
    if sum_x == 0 and sum_y == 0 and sum_z == 0:
        print("YES")
    else:
        print("NO")

solve()