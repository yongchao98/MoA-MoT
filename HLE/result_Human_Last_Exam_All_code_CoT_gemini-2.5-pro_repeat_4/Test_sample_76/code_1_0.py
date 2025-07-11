import math

def solve():
    """
    Finds the number of integers x for which x^3 - 16x^2 - 72x + 1056 is a perfect square.
    """
    solutions = []
    # We search in a range. For x < -8, the polynomial is negative.
    # For large x, the polynomial can be shown to lie between consecutive squares.
    # A search up to 10000 is more than sufficient.
    for x in range(-100, 10000):
        val = x**3 - 16*x**2 - 72*x + 1056
        if val >= 0:
            sqrt_val = math.isqrt(val)
            if sqrt_val * sqrt_val == val:
                solutions.append(x)
                print(f"For x = {x}, the expression is a perfect square: {x}^3 - 16*{x}^2 - 72*{x} + 1056 = {val} = {sqrt_val}^2")

    print(f"\nThe number of integer values of x for which the expression is a perfect square is: {len(solutions)}")

solve()
<<<3>>>