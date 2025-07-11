def check_equilibrium():
    """
    Reads n vectors of forces and determines if a body is in equilibrium.
    A body is in equilibrium if the vector sum of all forces is zero.
    """
    try:
        # Read the number of force vectors.
        n = int(input())
    except (ValueError, EOFError):
        # Handle invalid or empty input for n.
        n = 0

    # Initialize sums of force components for each axis.
    sum_x, sum_y, sum_z = 0, 0, 0

    # Loop n times to read each vector and accumulate the sums.
    for _ in range(n):
        try:
            # Read the three integer coordinates of a force vector.
            xi, yi, zi = map(int, input().split())
            
            # Add the components to the running totals.
            sum_x += xi
            sum_y += yi
            sum_z += zi
        except (ValueError, EOFError):
            # If a line is malformed or input ends prematurely,
            # the sums will likely be non-zero, leading to a "NO" result.
            # We can simply break or continue with potentially incorrect sums.
            # For this problem, we assume valid input as per the constraints.
            pass

    # Check if the body is in equilibrium.
    # This is true if the sum of forces along each axis is zero.
    if sum_x == 0 and sum_y == 0 and sum_z == 0:
        print("YES")
    else:
        print("NO")

check_equilibrium()