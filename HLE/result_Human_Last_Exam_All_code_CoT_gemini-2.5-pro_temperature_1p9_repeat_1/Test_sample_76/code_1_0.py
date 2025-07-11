import math

def solve():
    """
    This function finds the number of integers x for which
    P(x) = x^3 - 16x^2 - 72x + 1056 is a perfect square.
    
    The search for x is limited to the ranges [-8, 7] and [17, 62],
    as explained by the mathematical analysis.
    """
    
    solutions = []
    
    # Define the ranges to check
    ranges_to_check = [range(-8, 8), range(17, 63)]
    
    print("Searching for integer solutions x where P(x) = x^3 - 16x^2 - 72x + 1056 is a perfect square...")
    print("="*80)
    
    for r in ranges_to_check:
        for x in r:
            # Calculate P(x)
            p_x = x**3 - 16*x**2 - 72*x + 1056
            
            # P(x) must be non-negative
            if p_x >= 0:
                # Check if it's a perfect square
                sqrt_px = math.isqrt(p_x)
                if sqrt_px * sqrt_px == p_x:
                    solutions.append(x)
                    # Output the equation for each solution found
                    print(f"Found a solution for x = {x}:")
                    print(f"({x})^3 - 16*({x})^2 - 72*({x}) + 1056 = {p_x} = {sqrt_px}^2")
                    print("-" * 20)

    print("="*80)
    print(f"Total number of integers found: {len(solutions)}")

solve()
