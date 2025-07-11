import sys

def solve():
    """
    Determines if a body is in equilibrium based on a set of force vectors.
    """
    try:
        # Read the number of vectors
        n = int(sys.stdin.readline())
        
        sum_x = 0
        sum_y = 0
        sum_z = 0
        
        # Loop n times to read all vectors
        for _ in range(n):
            # Read a line and split it into components
            line = sys.stdin.readline()
            parts = line.split()
            
            # Convert parts to integers and add to sums
            xi = int(parts[0])
            yi = int(parts[1])
            zi = int(parts[2])
            
            sum_x += xi
            sum_y += yi
            sum_z += zi
            
        # Check for equilibrium condition
        if sum_x == 0 and sum_y == 0 and sum_z == 0:
            print("YES")
        else:
            print("NO")
            
    except (IOError, ValueError) as e:
        # Handle potential errors with input reading or conversion
        sys.stderr.write(f"An error occurred: {e}\n")

solve()