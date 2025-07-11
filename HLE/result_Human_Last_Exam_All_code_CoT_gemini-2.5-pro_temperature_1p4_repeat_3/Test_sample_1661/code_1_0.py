import sys

def solve_equilibrium_forces():
    """
    This function reads a series of force vectors and determines if a body
    subjected to these forces is in equilibrium.
    """
    try:
        # Read the number of force vectors from the first line of input.
        n = int(sys.stdin.readline())

        # Initialize variables to store the sum of forces for each axis.
        sum_x = 0
        sum_y = 0
        sum_z = 0

        # Loop 'n' times to read each force vector.
        for _ in range(n):
            # Read a line, split it by spaces, and convert the parts to integers.
            x, y, z = map(int, sys.stdin.readline().split())
            
            # Add the vector components to their respective sums.
            sum_x += x
            sum_y += y
            sum_z += z

        # A body is in equilibrium if the net force along each axis is zero.
        # This is checked by seeing if all three component sums are zero.
        if sum_x == 0 and sum_y == 0 and sum_z == 0:
            print("YES")
        else:
            print("NO")

    except (ValueError, IndexError):
        # This block catches potential errors, such as empty input or lines that
        # do not contain three valid integers. In such cases, the input is
        # malformed, so we can assume the body is not in equilibrium.
        print("NO")

# Execute the main function to solve the problem.
solve_equilibrium_forces()