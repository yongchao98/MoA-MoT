import sys

def check_equilibrium():
    """
    Reads n force vectors from standard input and determines if their
    sum is zero. A body is in equilibrium if the net force is zero.
    """
    try:
        # Read the number of vectors.
        n = int(sys.stdin.readline())
        
        # Initialize the components of the resultant force vector.
        sum_x, sum_y, sum_z = 0, 0, 0
        
        # Loop n times to read each vector and add it to the sum.
        for _ in range(n):
            # Read the three integer coordinates from a line.
            xi, yi, zi = map(int, sys.stdin.readline().split())
            
            # Add the components to the total sums.
            sum_x += xi
            sum_y += yi
            sum_z += zi
            
        # A body is in equilibrium if the sum of forces along each axis is zero.
        if sum_x == 0 and sum_y == 0 and sum_z == 0:
            print("YES")
        else:
            print("NO")

    except (ValueError, IndexError):
        # This handles cases of malformed input, though the problem
        # statement implies input will always be valid.
        print("NO")

check_equilibrium()