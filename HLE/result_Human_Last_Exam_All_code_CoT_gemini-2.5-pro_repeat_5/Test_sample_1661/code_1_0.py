import sys

def solve():
    """
    Determines if a body is in equilibrium based on a series of force vectors.
    """
    try:
        # Read the number of force vectors.
        n = int(sys.stdin.readline())
    except (ValueError, IndexError):
        # Handle cases with invalid or no input for n.
        print("NO")
        return

    # A body cannot be in equilibrium with 0 forces unless specified.
    # The problem states 1 <= n <= 100, so n will be at least 1.
    
    # Initialize sums for each coordinate.
    sum_x = 0
    sum_y = 0
    sum_z = 0

    # Loop n times to read each vector and add it to the sum.
    for _ in range(n):
        try:
            # Read a line, split it into parts, and convert to integers.
            x, y, z = map(int, sys.stdin.readline().split())
            
            # Add the components to the total force sums.
            sum_x += x
            sum_y += y
            sum_z += z
        except (ValueError, IndexError):
            # This handles malformed lines, though the problem description
            # implies input will always be valid.
            print("NO")
            return

    # A body is in equilibrium if the sum of force vectors is the zero vector.
    # This means the sum of components on each axis must be zero.
    if sum_x == 0 and sum_y == 0 and sum_z == 0:
        print("YES")
    else:
        print("NO")

solve()