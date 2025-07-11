import sys

def check_equilibrium():
    """
    Reads a series of 3D force vectors from standard input and determines
    if their sum is the zero vector, indicating the body is in equilibrium.
    """
    try:
        # Read the number of force vectors.
        n = int(sys.stdin.readline())
        
        # Ensure n is within the specified constraints.
        if not (1 <= n <= 100):
            # This case might not be strictly necessary based on problem statement,
            # but it's good practice. An invalid n would likely lead to "NO" anyway.
            print("NO")
            return

        sum_x, sum_y, sum_z = 0, 0, 0

        # Loop n times to read each vector and add it to the total sum.
        for _ in range(n):
            # Read the line and split it into components.
            line = sys.stdin.readline()
            parts = line.split()
            
            # Add components to their respective sums.
            sum_x += int(parts[0])
            sum_y += int(parts[1])
            sum_z += int(parts[2])

        # A body is in equilibrium if the net force on it is zero.
        # This means the sum of forces along each axis must be zero.
        if sum_x == 0 and sum_y == 0 and sum_z == 0:
            print("YES")
        else:
            print("NO")

    except (ValueError, IndexError):
        # If input is malformed (e.g., not enough numbers), it cannot be in equilibrium.
        print("NO")

check_equilibrium()