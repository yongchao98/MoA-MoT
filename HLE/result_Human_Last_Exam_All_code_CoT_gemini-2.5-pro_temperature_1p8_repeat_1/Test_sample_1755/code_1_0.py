import sys

def solve_equilibrium():
    """
    Reads force vectors from stdin, determines if the net force is zero,
    and prints the result.
    """
    try:
        # Read the number of vectors, n.
        line = sys.stdin.readline()
        if not line.strip():
            print("NO")
            return
        n = int(line)
        
        # Check constraints on n
        if not (1 <= n <= 100):
            # This case is not specified by the problem but is good practice.
            # Assuming valid n as per the problem description.
            pass

    except (ValueError, IndexError):
        # Invalid input for n
        print("NO")
        return

    # Initialize sums for x, y, z components.
    sum_x, sum_y, sum_z = 0, 0, 0
    
    # Loop n times to read each vector.
    for _ in range(n):
        try:
            line = sys.stdin.readline()
            if not line.strip():
                # Reached end of input unexpectedly.
                sum_x, sum_y, sum_z = 1, 1, 1 # Force a "NO" result
                break
            
            # Parse the three integer coordinates from the line.
            parts = list(map(int, line.split()))
            if len(parts) != 3:
                sum_x, sum_y, sum_z = 1, 1, 1 # Force a "NO" result
                break
                
            xi, yi, zi = parts
            
            # Add components to their respective sums.
            sum_x += xi
            sum_y += yi
            sum_z += zi
        except (ValueError, IndexError):
            # Line contains non-integer data or wrong number of parts.
            sum_x, sum_y, sum_z = 1, 1, 1 # Force a "NO" result
            break

    # The "final equation" involves checking if the sums are zero.
    # The C code uses the expression x || y || z.
    # In Python, this is equivalent to `sum_x != 0 or sum_y != 0 or sum_z != 0`.
    # Let's print the values used in this check.
    print(f"Final sums used in check: x={sum_x}, y={sum_y}, z={sum_z}")
    
    # Check the equilibrium condition and print the result.
    if sum_x == 0 and sum_y == 0 and sum_z == 0:
        print("Equation result: 0 || 0 || 0 => false")
        print("YES")
    else:
        # For non-programmers, x || y || z is tricky. Let's show the values.
        print(f"Equation result: {sum_x} || {sum_y} || {sum_z} => true")
        print("NO")

solve_equilibrium()