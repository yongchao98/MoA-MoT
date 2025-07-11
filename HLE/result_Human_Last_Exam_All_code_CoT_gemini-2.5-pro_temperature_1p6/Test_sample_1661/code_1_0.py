import sys

def solve():
    """
    Determines if a body is in equilibrium by summing force vectors.
    """
    try:
        # Read the number of force vectors
        n = int(sys.stdin.readline())
        
        # Initialize sums for each coordinate to zero
        sum_x = 0
        sum_y = 0
        sum_z = 0
        
        # Loop n times to read each vector and add it to the sums
        for _ in range(n):
            line = sys.stdin.readline()
            # In case of empty lines, skip
            if not line.strip():
                continue
            xi, yi, zi = map(int, line.split())
            sum_x += xi
            sum_y += yi
            sum_z += zi
            
        # As per instructions, output the numbers from the final equilibrium equation.
        # The equation is: ΣFx = 0, ΣFy = 0, ΣFz = 0. The numbers are the sums.
        print(f"Sum of forces: Fx={sum_x}, Fy={sum_y}, Fz={sum_z}")

        # Check for equilibrium condition
        if sum_x == 0 and sum_y == 0 and sum_z == 0:
            print("YES")
        else:
            print("NO")
            
    except (IOError, ValueError) as e:
        print(f"An error occurred with input: {e}", file=sys.stderr)

solve()