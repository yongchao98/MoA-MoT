import math

def is_perfect_square(n):
    """Checks if a number is a perfect square."""
    if n < 0:
        return False, -1
    if n == 0:
        return True, 0
    x = int(math.sqrt(n))
    return x * x == n, x

def p(x):
    """Calculates the value of the polynomial."""
    return x**3 - 16*x**2 - 72*x + 1056

def find_solutions():
    """Finds and prints integer solutions."""
    solutions = set()

    # Case 1: Search in the interval [-8, 16]
    for x in range(-8, 17):
        val = p(x)
        is_sq, y = is_perfect_square(val)
        if is_sq:
            solutions.add(x)

    # Case 2: x = k^2 for k up to 68
    for k in range(1, 69):
        x = k*k
        val = p(x)
        is_sq, y = is_perfect_square(val)
        if is_sq:
            solutions.add(x)

    # Case 3: x = m^2 + 16 for m up to 35
    for m in range(1, 36):
        x = m*m + 16
        val = p(x)
        is_sq, y = is_perfect_square(val)
        if is_sq:
            solutions.add(x)

    print("Integers x for which the quantity is a perfect square:")
    sorted_solutions = sorted(list(solutions))
    for x in sorted_solutions:
        val = p(x)
        is_sq, y = is_perfect_square(val)
        print(f"x = {x}, P(x) = {val} = {y}^2")
    
    print(f"\nTotal number of integers found: {len(sorted_solutions)}")

if __name__ == '__main__':
    find_solutions()