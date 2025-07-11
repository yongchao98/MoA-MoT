import sys

def solve():
    """
    Reads n vectors and determines if their sum is the zero vector.
    """
    try:
        # Read the number of vectors, n.
        n = int(sys.stdin.readline())
        
        # Initialize sums for each coordinate.
        sum_x, sum_y, sum_z = 0, 0, 0
        
        # Loop n times to read each vector.
        for _ in range(n):
            line = sys.stdin.readline()
            if not line:
                break 
            
            parts = list(map(int, line.split()))
            if len(parts) == 3:
                x, y, z = parts
                sum_x += x
                sum_y += y
                sum_z += z

        # Check if the body is in equilibrium.
        # This is true if and only if the sum of vectors is the zero vector.
        if sum_x == 0 and sum_y == 0 and sum_z == 0:
            print("YES")
        else:
            # We also print the final sum vector for clarity
            # print(f"Final vector sum: <{sum_x}, {sum_y}, {sum_z}>")
            print("NO")

    except (ValueError, IndexError):
        # Handle cases with bad input format.
        print("NO")

solve()