import math

def find_perfect_square_solutions():
    """
    This function finds all integers x for which the expression
    x^3 - 16x^2 - 72x + 1056 is a perfect square.
    """
    solutions = []
    # Based on mathematical analysis, we only need to check the integer range [-8, 62].
    for x in range(-8, 63):
        val = x**3 - 16*x**2 - 72*x + 1056
        if val >= 0:
            sqrt_val = math.isqrt(val)
            if sqrt_val * sqrt_val == val:
                solutions.append((x, sqrt_val))

    print(f"Found {len(solutions)} integer(s) for x where the expression is a perfect square.\n")

    if solutions:
        print("The solutions are:")
        for x, y in solutions:
            p1 = x**3
            p2 = -16 * (x**2)
            p3 = -72 * x
            p4 = 1056
            
            # Formatting the output to show each number in the equation.
            print(f"For x = {x}:")
            print(f"({x})^3 - 16*({x})^2 - 72*({x}) + 1056 = {p1} + ({p2}) + ({p3}) + {p4} = {val} = {y}^2")
            print("-" * 20)

find_perfect_square_solutions()
