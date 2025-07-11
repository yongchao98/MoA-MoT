import math

def is_perfect_square(n):
    """Checks if a number is a perfect square."""
    if n < 0:
        return False
    if n == 0:
        return True
    x = int(math.sqrt(n))
    return x * x == n

def p(x):
    """Calculates the value of the polynomial for a given x."""
    return x**3 - 16*x**2 - 72*x + 1056

def find_solutions():
    """Finds and prints all integer solutions in the determined range."""
    solutions = []
    # Based on analysis, we only need to check x in [-38, 62].
    for x in range(-38, 63):
        val = p(x)
        if is_perfect_square(val):
            solutions.append(x)

    print("The integer values of x for which the expression is a perfect square are:")
    if not solutions:
        print("None found in the search range.")
    else:
        for x in solutions:
            y_squared = p(x)
            y = int(math.sqrt(y_squared))
            term1 = x**3
            term2 = -16 * x**2
            term3 = -72 * x
            term4 = 1056
            # The prompt asks to output each number in the final equation.
            print(f"For x = {x}:")
            print(f"({x})^3 - 16*({x})^2 - 72*({x}) + 1056 = {term1} + ({term2}) + ({term3}) + {term4} = {y_squared} = {y}^2")
    
    print(f"\nThe number of such integers is {len(solutions)}.")

if __name__ == '__main__':
    find_solutions()