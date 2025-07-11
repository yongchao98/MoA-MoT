import sys

def check_equilibrium():
    """
    Reads n 3D force vectors and determines if their sum is zero.
    """
    try:
        # Read the number of vectors
        n_str = sys.stdin.readline()
        if not n_str:
            return
        n = int(n_str)
        
        sum_x = 0
        sum_y = 0
        sum_z = 0
        
        # Loop n times to read each vector and add it to the sum
        for _ in range(n):
            line = sys.stdin.readline()
            if not line:
                break
            
            # Split the line into parts and convert to integers
            parts = line.split()
            x = int(parts[0])
            y = int(parts[1])
            z = int(parts[2])
            
            # Add to the component-wise sums
            sum_x += x
            sum_y += y
            sum_z += z
            
        # Check for equilibrium condition
        if sum_x == 0 and sum_y == 0 and sum_z == 0:
            print("YES")
        else:
            print("NO")

    except (ValueError, IndexError):
        # Handle potential errors with input format
        print("NO")

check_equilibrium()