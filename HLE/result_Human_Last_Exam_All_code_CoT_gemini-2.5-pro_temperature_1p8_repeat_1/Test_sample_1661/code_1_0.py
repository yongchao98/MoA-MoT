import sys

def solve():
    """
    Reads force vectors and determines if the body is in equilibrium.
    """
    try:
        # Read the number of vectors
        n = int(sys.stdin.readline())
    except (ValueError, IndexError):
        # Handle cases where input for n is invalid or missing
        n = 0

    sum_x = 0
    sum_y = 0
    sum_z = 0

    # Loop n times to read each vector
    for _ in range(n):
        try:
            # Read a line, split it by spaces, and convert components to integers
            line = sys.stdin.readline().strip()
            if not line:
                continue
            xi, yi, zi = map(int, line.split())
            
            # Add the components to the total sums
            sum_x += xi
            sum_y += yi
            sum_z += zi
        except (ValueError, IndexError):
            # Assuming valid input format as per problem description
            pass
            
    # Per the instructions, "output each number in the final equation!"
    # The final equation is that the resultant force vector <sum_x, sum_y, sum_z> is zero.
    # The numbers are the components of that resultant vector.
    # Note: This line is not part of the original problem's output requirements,
    # but is included based on the assistant's instructions.
    # Let's format it in a simple way that just shows the numbers.
    print(f"{sum_x} {sum_y} {sum_z}")


    # A body is in equilibrium if the net force is zero in all directions
    if sum_x == 0 and sum_y == 0 and sum_z == 0:
        print("YES")
    else:
        print("NO")

solve()