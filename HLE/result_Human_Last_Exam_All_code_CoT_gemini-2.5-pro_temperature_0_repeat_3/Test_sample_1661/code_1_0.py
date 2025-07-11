import sys

def solve():
    """
    Determines if a body is in equilibrium based on a set of force vectors.
    """
    try:
        # Read the number of force vectors.
        n = int(sys.stdin.readline())
        
        # Initialize sums for each coordinate to zero.
        sum_x = 0
        sum_y = 0
        sum_z = 0
        
        # Loop n times to read each vector and update the sums.
        for _ in range(n):
            line = sys.stdin.readline()
            # Read the three integer coordinates from the line.
            x, y, z = map(int, line.split())
            sum_x += x
            sum_y += y
            sum_z += z
            
        # Check if the body is in equilibrium.
        # This occurs if the sum of vectors on each axis is 0.
        if sum_x == 0 and sum_y == 0 and sum_z == 0:
            print("YES")
        else:
            print("NO")
            
    except (IOError, ValueError) as e:
        # Handle potential input errors.
        print(f"An error occurred: {e}", file=sys.stderr)

solve()