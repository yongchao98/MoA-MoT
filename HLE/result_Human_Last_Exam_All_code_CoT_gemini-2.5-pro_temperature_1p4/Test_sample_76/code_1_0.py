import math

def is_perfect_square(n):
    """Checks if a number is a perfect square."""
    if n < 0:
        return False
    if n == 0:
        return True
    # Use integer square root to avoid floating point inaccuracies
    s = int(math.sqrt(n))
    return s * s == n

def find_integer_solutions():
    """
    Finds the integer values of x for which the expression
    x^3 - 16x^2 - 72x + 1056 is a perfect square.
    """
    solutions = set()
    
    # Based on analysis, we check integers x in [-5, 4] and [17, 62]
    # In these ranges, P(x) can be non-negative.
    candidate_ranges = list(range(-5, 5)) + list(range(17, 63))
    
    for x in candidate_ranges:
        val = x**3 - 16*x**2 - 72*x + 1056
        if is_perfect_square(val):
            solutions.add(x)
            
    # For x >= 63, our analysis shows x must be a perfect square, k^2.
    # For k >= 69, analysis shows there are no solutions.
    # So we check k from ceil(sqrt(63))=8 to 68.
    for k in range(8, 69):
        x = k*k
        val = x**3 - 16*x**2 - 72*x + 1056
        if is_perfect_square(val):
            solutions.add(x)
            
    # Sort the found solutions
    solution_list = sorted(list(solutions))
    
    # Print the results
    print(f"Found {len(solution_list)} integer solution(s) for x.")
    if solution_list:
        print("The solutions are:")
        for x in solution_list:
            y_squared = x**3 - 16*x**2 - 72*x + 1056
            y = int(math.sqrt(y_squared))
            
            # Deconstruct the calculation for clarity
            term1 = x**3
            term2 = -16 * (x**2)
            term3 = -72 * x
            term4 = 1056
            
            print(f"For x = {x}:")
            print(f"({x})^3 - 16*({x})^2 - 72*({x}) + 1056 = {term1} + ({term2}) + ({term3}) + {term4} = {y_squared}, which is {y}^2.")

if __name__ == '__main__':
    find_integer_solutions()
