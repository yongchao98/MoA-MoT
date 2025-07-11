import sys

def solve():
    """
    Determines if a body is in equilibrium based on n force vectors.
    """
    try:
        # Read the number of vectors, n.
        n = int(sys.stdin.readline())
        
        # Initialize sums of vector components.
        sum_x = 0
        sum_y = 0
        sum_z = 0
        
        # Loop n times to read all vectors.
        for _ in range(n):
            line = sys.stdin.readline()
            if not line:
                break
            # Read the three components, convert them to int, and add to sums.
            try:
                xi, yi, zi = map(int, line.split())
                sum_x += xi
                sum_y += yi
                sum_z += zi
            except (ValueError, IndexError):
                # Handle potential empty or malformed lines
                continue
        
        # Check for equilibrium condition.
        if sum_x == 0 and sum_y == 0 and sum_z == 0:
            print("YES")
        else:
            print("NO")

    except (ValueError, IndexError):
        # Handle cases where the first line (n) is missing or invalid.
        print("NO")

solve()
