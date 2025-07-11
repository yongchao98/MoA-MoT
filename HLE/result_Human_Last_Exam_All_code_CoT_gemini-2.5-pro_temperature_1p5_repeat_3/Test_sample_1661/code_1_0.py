import sys

def check_equilibrium():
    """
    Reads a set of 3D force vectors from standard input and determines
    if they result in equilibrium (i.e., the net force is a zero vector).

    Input format:
    The first line contains an integer n, the number of vectors.
    The next n lines each contain three space-separated integers (x, y, z)
    representing a force vector.
    """
    try:
        # Read the number of force vectors.
        # The input is guaranteed to be a positive integer.
        num_vectors = int(sys.stdin.readline())

        # Initialize sums of vector components to zero.
        # The body is in equilibrium if the sum of all force vectors is the zero vector.
        # This means the sum of the x, y, and z components must each be zero.
        sum_x, sum_y, sum_z = 0, 0, 0

        # Loop n times to read each vector and add its components to the sums.
        for _ in range(num_vectors):
            line = sys.stdin.readline()
            # Split the line into components and convert them to integers.
            x, y, z = map(int, line.split())
            
            # Add the components to their respective sums.
            sum_x += x
            sum_y += y
            sum_z += z

        # After summing all vectors, check if the resultant vector is zero.
        if sum_x == 0 and sum_y == 0 and sum_z == 0:
            print("YES")
        else:
            print("NO")

    except (ValueError, IndexError):
        # This handles cases of malformed input lines.
        sys.stderr.write("Error: Invalid input format.\n")

# Execute the main function.
check_equilibrium()