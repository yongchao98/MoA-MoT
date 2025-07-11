import math

def solve():
    """
    Calculates the minimum value of 1000m + n based on the derived polynomial degree.
    """
    # Step 1: Determine the degree of the polynomial P.
    # Based on our analysis, the condition is the union of two hypersurfaces:
    # 1. A, B, C, D are coplanar: degree d1 = 3.
    # 2. A, B, C, D, X are cospherical: degree d2 = 5.
    # The total degree of the polynomial P is the sum of the degrees.
    d = 3 + 5
    print(f"The total degree of the polynomial is d = {d}")

    min_val = float('inf')
    best_m = -1
    best_n = -1

    # Step 2: Find the minimum of 1000m + n subject to n*m >= d.
    # We iterate through possible values of m from 1 to d.
    for m in range(1, d + 1):
        # For a given m, the smallest integer n must satisfy n*m >= d.
        n = math.ceil(d / m)
        
        current_val = 1000 * m + n
        
        if current_val < min_val:
            min_val = current_val
            best_m = m
            best_n = n

    # Step 3: Output the results.
    # The prompt asks to output each number in the final equation.
    print(f"To minimize 1000m + n for a polynomial of degree {d}:")
    print(f"The optimal value for m is {best_m}.")
    print(f"The optimal value for n is {best_n}.")
    print(f"The final equation is: 1000 * {best_m} + {best_n} = {min_val}")
    print(f"The minimum value is {min_val}.")

solve()