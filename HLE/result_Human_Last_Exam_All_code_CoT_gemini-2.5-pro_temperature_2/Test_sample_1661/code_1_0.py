import sys

def check_equilibrium():
    """
    Reads n force vectors and determines if their sum is the zero vector.
    A body is in equilibrium if the net force acting on it is zero.
    """
    try:
        # 1. Read the number of vectors
        n = int(sys.stdin.readline())

        # 2. Initialize sum of vector components to zero
        sum_x, sum_y, sum_z = 0, 0, 0

        # 3. Loop n times to read all vectors and sum their components
        for _ in range(n):
            line = sys.stdin.readline()
            components = list(map(int, line.split()))
            if len(components) == 3:
                sum_x += components[0]
                sum_y += components[1]
                sum_z += components[2]

        # Note on the "final equation": The equilibrium condition is the equation
        # F_net = (sum_x)i + (sum_y)j + (sum_z)k = 0.
        # The C code uses the numbers `x[0]`, `x[1]`, `x[2]` (our `sum_x`, `sum_y`, `sum_z`)
        # in the final logical check. Since the problem statement only asks for a "YES" or "NO"
        # output, we will stick to that format. The numbers are implicitly used in the final check.

        # 5. Check if all sums are zero and print the result
        if sum_x == 0 and sum_y == 0 and sum_z == 0:
            print("YES")
        else:
            print("NO")

    except (ValueError, IndexError):
        # Handle cases with malformed input
        print("Invalid input. Please check the format.")

# Execute the function
check_equilibrium()
