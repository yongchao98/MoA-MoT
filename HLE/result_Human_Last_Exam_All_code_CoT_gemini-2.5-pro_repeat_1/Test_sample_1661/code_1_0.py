import sys

def solve():
    """
    Determines if a body is in equilibrium based on force vectors.
    """
    try:
        # Read the number of force vectors.
        n = int(sys.stdin.readline())
    except (ValueError, IndexError):
        # Handle cases with no input or invalid input for n.
        print("NO")
        return

    # Initialize the sum of forces for each axis.
    sum_x = 0
    sum_y = 0
    sum_z = 0

    # Read n vectors and calculate the sum of components.
    for _ in range(n):
        try:
            line = sys.stdin.readline()
            # Parse the three integer components from the line.
            xi, yi, zi = map(int, line.split())
            
            # Add the components to the total sums.
            sum_x += xi
            sum_y += yi
            sum_z += zi
        except (ValueError, IndexError):
            # If a line is malformed or empty, we can treat its contribution as zero
            # or handle as an error. For this problem, we'll continue.
            continue
            
    # A body is in equilibrium if the net force vector is (0, 0, 0).
    if sum_x == 0 and sum_y == 0 and sum_z == 0:
        print("YES")
    else:
        print("NO")

solve()