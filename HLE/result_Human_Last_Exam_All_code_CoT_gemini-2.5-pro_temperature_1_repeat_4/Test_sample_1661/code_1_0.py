import sys

def check_equilibrium():
    """
    Reads n force vectors from standard input and determines if their
    sum is zero, indicating the body is in equilibrium.
    """
    try:
        # Read the number of vectors.
        n = int(sys.stdin.readline())

        # Initialize the sum of force components to zero.
        sum_x = 0
        sum_y = 0
        sum_z = 0

        # Read each vector and add its components to the sums.
        for _ in range(n):
            components = list(map(int, sys.stdin.readline().split()))
            sum_x += components[0]
            sum_y += components[1]
            sum_z += components[2]

        # A body is in equilibrium if the net force is the zero vector.
        # This means the sum of components along each axis must be zero.
        if sum_x == 0 and sum_y == 0 and sum_z == 0:
            print("YES")
        else:
            print("NO")

    except (ValueError, IndexError):
        # This handles cases with malformed input, though the problem
        # description implies the input will always be valid.
        print("NO")

check_equilibrium()