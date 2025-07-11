import math

def solve():
    """
    This function finds the number of integers x for which the expression
    x^3 - 16x^2 - 72x + 1056 is a perfect square.
    """
    solutions = []
    
    # The search range is determined by analyzing the polynomial.
    # P(x) >= 0 for x in [-8, 7] and x >= 17.
    # A further analysis shows no solutions exist for x > 110.
    ranges_to_check = [range(-8, 8), range(17, 111)]

    for r in ranges_to_check:
        for x in r:
            # Calculate the value of the polynomial
            p_x = x**3 - 16*x**2 - 72*x + 1056

            # A perfect square must be non-negative
            if p_x >= 0:
                # Check if p_x is a perfect square
                sqrt_px = math.isqrt(p_x)
                if sqrt_px * sqrt_px == p_x:
                    solutions.append((x, sqrt_px))

    print("The integer solutions (x, y) for y^2 = x^3 - 16x^2 - 72x + 1056 are:")
    for x, y in solutions:
        print(f"For x = {x}, the expression is {y}^2 = {p_x}")

    print("\nThe equations for the found integer solutions are:")
    for x, y in solutions:
        print(f"{x}^3 - 16*({x})^2 - 72*({x}) + 1056 = {y**2}")

    print(f"\nThe total number of such integers x is: {len(solutions)}")

if __name__ == "__main__":
    solve()
<<<3>>>