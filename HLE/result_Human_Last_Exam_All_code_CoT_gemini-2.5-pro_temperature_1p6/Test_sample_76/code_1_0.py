import math

def find_integer_solutions():
    """
    Finds the number of integers x for which x^3 - 16x^2 - 72x + 1056
    is a perfect square.
    """
    solutions = []
    
    # Define the polynomial P(x)
    def p(x):
        return x**3 - 16*x**2 - 72*x + 1056

    # Function to check if a number is a perfect square
    def is_perfect_square(n):
        if n < 0:
            return False
        sqrt_n = int(math.sqrt(n))
        return sqrt_n * sqrt_n == n

    # First range of interest: x in [-6, 4]
    for x in range(-6, 5):
        val = p(x)
        if is_perfect_square(val):
            solutions.append(x)

    # Second range of interest: x >= 18
    # We search up to a large bound, e.g., 30000, as theory suggests
    # solutions are finite and "small".
    for x in range(18, 30001):
        val = p(x)
        if is_perfect_square(val):
            solutions.append(x)
            
    print("The integer solutions for x are:")
    if not solutions:
        print("None")
    else:
        for sol in solutions:
            y = int(math.sqrt(p(sol)))
            print(f"x = {sol}, gives {sol}^3 - 16*{sol}^2 - 72*{sol} + 1056 = {y}^2")

    print(f"\nThe number of integers x is {len(solutions)}.")

if __name__ == '__main__':
    find_integer_solutions()