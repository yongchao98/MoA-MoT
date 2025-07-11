import sys

def solve():
    """
    Determines if a body is in equilibrium based on n force vectors.
    """
    try:
        # Read the number of vectors from the first line of input.
        n = int(sys.stdin.readline())
    except (ValueError, IndexError):
        # Handle case where input is empty or not a valid integer.
        print("NO")
        return

    # Initialize variables to store the sum of components for each axis.
    sum_x = 0
    sum_y = 0
    sum_z = 0

    # Loop n times to read each force vector.
    for _ in range(n):
        try:
            # Read a line, split it by spaces, and convert parts to integers.
            x, y, z = map(int, sys.stdin.readline().split())
            
            # Add the components to their respective sums.
            sum_x += x
            sum_y += y
            sum_z += z
        except (ValueError, IndexError):
            # If a line is malformed, we can assume the net force is not zero
            # or handle it as per specific requirements. Here we'll print NO
            # and exit for robustness.
            print("NO")
            return

    # A body is in equilibrium (idle) if the sum of forces along each axis is zero.
    # The condition for equilibrium is F_net = 0, which means ΣFx = 0, ΣFy = 0, and ΣFz = 0.
    if sum_x == 0 and sum_y == 0 and sum_z == 0:
        print("YES")
    else:
        print("NO")

solve()