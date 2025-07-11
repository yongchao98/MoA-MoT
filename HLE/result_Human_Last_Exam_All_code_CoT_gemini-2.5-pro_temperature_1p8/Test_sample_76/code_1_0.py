import math

def is_perfect_square(n):
    """Checks if a number is a perfect square."""
    if n < 0:
        return False
    if n == 0:
        return True
    root = int(math.sqrt(n))
    return root * root == n

def solve_problem():
    """
    Finds the number of integers x for which x^3 - 16x^2 - 72x + 1056 is a perfect square.
    """
    solutions = []
    
    # We check the ranges [-6, 7] and [17, 144] based on mathematical analysis.
    ranges_to_check = list(range(-6, 8)) + list(range(17, 145))
    
    for x in ranges_to_check:
        # Define the expression
        val = x**3 - 16*x**2 - 72*x + 1056
        if is_perfect_square(val):
            solutions.append(x)
            root = int(math.sqrt(val))
            # Outputting each number in the final equation
            # For negative numbers, we use parenthesis to make it clear.
            if x < 0:
                x_str = f"({x})"
            else:
                x_str = str(x)
                
            print(f"For x = {x}, the expression is {x_str}^3 - 16*{x_str}^2 - 72*{x_str} + 1056 = {val} = {root}^2")

    print(f"\nIn total, there are {len(solutions)} integer values of x for which the expression is a perfect square.")

solve_problem()