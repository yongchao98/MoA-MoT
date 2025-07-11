def check_equilibrium():
    """
    Reads n 3D force vectors and determines if their sum is zero.
    """
    try:
        # Read the number of vectors.
        n = int(input())
    except (ValueError, EOFError):
        # Handle invalid or no input for n.
        # An object with no forces is in equilibrium.
        # But based on constraints 1 <= n <= 100, this case shouldn't occur.
        # We'll print NO for safety.
        print("NO")
        return

    # Initialize sums for each component of the net force vector.
    sum_x = 0
    sum_y = 0
    sum_z = 0

    # Read each of the n vectors and update the sums.
    for _ in range(n):
        try:
            # Read the three components of a force vector.
            x, y, z = map(int, input().split())
            
            # Add the components to the running totals.
            sum_x += x
            sum_y += y
            sum_z += z
        except (ValueError, EOFError):
            # If a line of input is malformed, we can assume equilibrium is broken.
            print("NO")
            return
            
    # A body is in equilibrium if the sum of forces along each axis is zero.
    if sum_x == 0 and sum_y == 0 and sum_z == 0:
        print("YES")
    else:
        print("NO")

check_equilibrium()
