import math

def is_perfect_square(n):
    """
    Checks if a non-negative integer is a perfect square.
    """
    if n < 0:
        return False
    if n == 0:
        return True
    x = int(math.sqrt(n))
    return x * x == n

def find_integer_solutions():
    """
    Finds and prints the integer values of x for which the expression
    x^3 - 16x^2 - 72x + 1056 is a perfect square.
    
    The problem is to find integer solutions to y^2 = x^3 - 16x^2 - 72x + 1056.
    Analysis shows that we only need to search for solutions in the integer ranges
    x in [-8, 7] and x in [17, 62].
    """
    solutions = []
    
    # Define the ranges to check based on prior analysis
    search_ranges = [range(-8, 8), range(17, 63)]
    
    for r in search_ranges:
        for x in r:
            # The polynomial is P(x) = x^3 - 16x^2 - 72x + 1056
            val = x**3 - 16*x**2 - 72*x + 1056
            
            if is_perfect_square(val):
                y = int(math.sqrt(val))
                solutions.append({'x': x, 'val': val, 'y': y})

    print(f"The number of integers x for which the quantity is a perfect square is {len(solutions)}.")
    
    if solutions:
        print("The solutions are:")
        for s in solutions:
            x = s['x']
            val = s['val']
            y = s['y']
            
            # Calculate each term to display in the equation
            term1 = x**3
            term2 = -16 * x**2
            term3 = -72 * x
            term4 = 1056
            
            print(f"For x = {x}:")
            print(f"({x})^3 + ({term2}) + ({term3}) + {term4} = {val} = {y}^2")

# Run the solver
find_integer_solutions()