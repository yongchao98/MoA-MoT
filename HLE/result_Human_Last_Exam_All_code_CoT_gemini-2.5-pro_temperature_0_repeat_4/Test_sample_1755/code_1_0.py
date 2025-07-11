import sys

def solve_equilibrium():
    """
    Reads force vectors from standard input and determines if a body is in equilibrium.
    """
    try:
        # Read the number of vectors, n.
        line = sys.stdin.readline()
        if not line.strip():
            return
        n = int(line)
    except (ValueError, IndexError):
        # Handle cases where n is not provided or is not a valid integer.
        print("NO")
        return

    # Initialize sums of vector components
    sum_x = 0
    sum_y = 0
    sum_z = 0

    # Loop n times to read each vector
    for _ in range(n):
        line = sys.stdin.readline()
        if not line.strip():
            # Stop if input ends prematurely
            break
        try:
            # Split the line, convert components to integers, and unpack
            xi, yi, zi = map(int, line.split())
            
            # Add components to the total sums
            sum_x += xi
            sum_y += yi
            sum_z += zi
        except (ValueError, IndexError):
            # If a line is malformed, it could be considered as contributing nothing,
            # but for this problem, we assume valid input as per the constraints.
            # We'll treat malformed lines as an error condition leading to "NO".
            sum_x = 1 # Set to a non-zero value to ensure the result is "NO"
            break

    # A body is in equilibrium if the sum of all force vectors is the zero vector.
    if sum_x == 0 and sum_y == 0 and sum_z == 0:
        print("YES")
    else:
        print("NO")

solve_equilibrium()