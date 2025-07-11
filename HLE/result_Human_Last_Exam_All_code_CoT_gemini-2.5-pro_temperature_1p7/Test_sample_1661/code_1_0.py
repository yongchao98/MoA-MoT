import sys

def solve_equilibrium():
    """
    Reads n force vectors and determines if their sum is zero.
    A body is in equilibrium if the net force acting on it is the zero vector.
    This means the sum of the x-components, y-components, and z-components
    of all force vectors must each be zero.
    """
    try:
        # Read the number of force vectors from the first line of input.
        num_vectors = int(sys.stdin.readline())

        # Initialize the sums of the vector components to zero.
        sum_x, sum_y, sum_z = 0, 0, 0

        # Loop n times to read each vector and update the sums.
        for _ in range(num_vectors):
            # Read a line, split it into parts, and convert them to integers.
            parts = sys.stdin.readline().split()
            x, y, z = map(int, parts)
            
            # Add the components of the current vector to the total sums.
            sum_x += x
            sum_y += y
            sum_z += z

        # Check if the body is in equilibrium.
        # This occurs if and only if the sum of components for each axis is zero.
        if sum_x == 0 and sum_y == 0 and sum_z == 0:
            print("YES")
        else:
            print("NO")

    except (ValueError, IndexError):
        # According to the problem description, the input will always be valid.
        # This handles potential errors in other scenarios.
        print("NO")

solve_equilibrium()