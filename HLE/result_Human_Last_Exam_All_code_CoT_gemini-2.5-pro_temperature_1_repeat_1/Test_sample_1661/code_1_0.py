import sys

def solve():
    """
    Determines if a body is in equilibrium based on the sum of force vectors.
    """
    try:
        # Read the number of force vectors
        n = int(sys.stdin.readline())
        
        # Initialize sums of force components
        sum_x = 0
        sum_y = 0
        sum_z = 0
        
        # Loop n times to read all vectors and sum their components
        for _ in range(n):
            line = sys.stdin.readline()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            
            xi = int(parts[0])
            yi = int(parts[1])
            zi = int(parts[2])
            
            sum_x += xi
            sum_y += yi
            sum_z += zi
            
        # Check if the body is in equilibrium
        if sum_x == 0 and sum_y == 0 and sum_z == 0:
            print("YES")
        else:
            print("NO")
            
    except (IOError, ValueError) as e:
        # Handle potential input errors
        print(f"An error occurred: {e}", file=sys.stderr)

solve()