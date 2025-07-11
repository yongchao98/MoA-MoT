import math

def solve():
    """
    This function calculates the number of integers in the specified range
    that can be expressed in the given form.
    """
    lower_bound = 10**18
    upper_bound = 10**18 + 10000
    
    solutions = set()

    # Case 1: y = 0, x = 10^6
    # n = 10^18 + 4z^3
    # We need 0 <= 4z^3 <= 10000, which means 0 <= z^3 <= 2500.
    z_max = int(math.pow(2500, 1/3))
    print("Solutions for n = 10^18 + 4*z^3:")
    for z in range(z_max + 1):
        # For z to be an integer, z must be non-negative.
        n = lower_bound + 4 * (z**3)
        if lower_bound <= n <= upper_bound:
            solutions.add(n)
            print(f"For z = {z}: n = 10^18 + 4*({z})^3 = 10^18 + {4 * z**3}")

    # Case 2: z = 0, x = 10^6
    # n = 10^18 + 2y^3
    # We need 0 <= 2y^3 <= 10000, which means 0 <= y^3 <= 5000.
    y_max = int(math.pow(5000, 1/3))
    print("\nSolutions for n = 10^18 + 2*y^3:")
    for y in range(y_max + 1):
        # For y to be an integer, y must be non-negative.
        n = lower_bound + 2 * (y**3)
        if lower_bound <= n <= upper_bound:
            solutions.add(n)
            print(f"For y = {y}: n = 10^18 + 2*({y})^3 = 10^18 + {2 * y**3}")

    print(f"\nTotal number of distinct integers is: {len(solutions)}")

solve()