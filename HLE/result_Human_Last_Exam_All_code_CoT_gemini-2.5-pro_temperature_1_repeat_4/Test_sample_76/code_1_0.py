import math

def is_perfect_square(n):
    """Checks if a number is a perfect square."""
    if n < 0:
        return False
    if n == 0:
        return True
    x = int(math.sqrt(n))
    return x * x == n

def find_solutions():
    """Finds and prints integer solutions for x."""
    solutions = []
    
    # Check the first interval [-8, 7]
    for x in range(-8, 8):
        val = x**3 - 16*x**2 - 72*x + 1056
        if is_perfect_square(val):
            solutions.append(x)
            y = int(math.sqrt(val))
            print(f"For x = {x}:")
            print(f"{x}^3 - 16*({x})^2 - 72*({x}) + 1056 = {val} = {y}^2")

    # Check the second interval [17, ...]
    # A search up to 5000 is sufficient based on bounding analysis.
    for x in range(17, 5001):
        val = x**3 - 16*x**2 - 72*x + 1056
        if is_perfect_square(val):
            solutions.append(x)
            y = int(math.sqrt(val))
            print(f"For x = {x}:")
            print(f"{x}^3 - 16*({x})^2 - 72*({x}) + 1056 = {val} = {y}^2")
    
    print(f"\nThe integer values of x are: {', '.join(map(str, sorted(solutions)))}.")
    print(f"The number of such integers is {len(solutions)}.")

if __name__ == '__main__':
    find_solutions()
