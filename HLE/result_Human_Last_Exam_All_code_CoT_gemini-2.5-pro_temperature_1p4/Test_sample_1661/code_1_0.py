import sys

def solve():
    """
    Determines if a body is in equilibrium by summing force vectors.
    """
    try:
        # Read the number of force vectors.
        n = int(sys.stdin.readline())
        
        # Initialize sums for each vector component.
        sum_x = 0
        sum_y = 0
        sum_z = 0
        
        # Loop n times to read each vector and add it to the sum.
        for _ in range(n):
            line = sys.stdin.readline()
            parts = line.split()
            if len(parts) == 3:
                fx = int(parts[0])
                fy = int(parts[1])
                fz = int(parts[2])
                sum_x += fx
                sum_y += fy
                sum_z += fz
        
        # Check if the net force is the zero vector.
        if sum_x == 0 and sum_y == 0 and sum_z == 0:
            print("YES")
        else:
            print("NO")
            
    except (ValueError, IndexError) as e:
        # Handle potential errors with input format.
        print(f"Error processing input: {e}", file=sys.stderr)

solve()