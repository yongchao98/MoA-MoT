import sys

def check_equilibrium():
    """
    Determines if a body is in equilibrium by summing force vectors.
    """
    try:
        # Read the number of force vectors from standard input.
        n = int(sys.stdin.readline())

        # According to the problem constraints, 1 <= n <= 100.
        # If n is out of this range, the behavior is undefined by the problem,
        # but we'll proceed.

        # Initialize the sums of the x, y, and z components of the forces.
        sum_x = 0
        sum_y = 0
        sum_z = 0

        # Loop n times to read each force vector and accumulate the components.
        for _ in range(n):
            # Read a line, split it into parts, and convert them to integers.
            parts = sys.stdin.readline().split()
            x = int(parts[0])
            y = int(parts[1])
            z = int(parts[2])

            # Add the components to their respective sums.
            sum_x += x
            sum_y += y
            sum_z += z

        # The final equation for equilibrium is that the resultant force vector is (0, 0, 0).
        # This means sum_x must be 0, sum_y must be 0, and sum_z must be 0.
        # The numbers in this final equation are the calculated sums.
        # The problem asks only for a "YES" or "NO" output, not the sums themselves.
        
        # Check if the body is in equilibrium.
        if sum_x == 0 and sum_y == 0 and sum_z == 0:
            print("YES")
        else:
            print("NO")

    except (ValueError, IndexError):
        # This handles cases of malformed input, such as non-integer values
        # or an incorrect number of components for a vector.
        # A robust program should handle this; for the specified problem,
        # we can assume valid input, but printing "NO" is a safe default.
        print("NO")

# Execute the main function.
check_equilibrium()