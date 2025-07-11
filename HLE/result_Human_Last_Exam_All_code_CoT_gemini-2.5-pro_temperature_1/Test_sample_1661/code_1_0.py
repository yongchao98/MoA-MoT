import sys

def solve_equilibrium():
    """
    Determines if a body is in equilibrium given n force vectors.
    A body is in equilibrium if the net force is zero, which means the
    sum of vector components in each dimension (x, y, z) is zero.
    """
    try:
        # Read the number of force vectors.
        n = int(sys.stdin.readline())
        
        if not (1 <= n <= 100):
            # Optional: Handle constraint violation, though problem statement implies valid input.
            return

        # Initialize sums for x, y, and z components.
        sum_x, sum_y, sum_z = 0, 0, 0

        # Read n vectors and accumulate the sums.
        for _ in range(n):
            line = sys.stdin.readline()
            parts = line.split()
            # It's good practice to ensure there are 3 parts.
            if len(parts) == 3:
                fx = int(parts[0])
                fy = int(parts[1])
                fz = int(parts[2])
                sum_x += fx
                sum_y += fy
                sum_z += fz
        
        # Check if the body is in equilibrium.
        if sum_x == 0 and sum_y == 0 and sum_z == 0:
            print("YES")
        else:
            print("NO")

    except (ValueError, IndexError):
        # Handle potential errors from input conversion or splitting.
        print("NO") # Assuming invalid input means it's not in a stable equilibrium.

if __name__ == "__main__":
    solve_equilibrium()