import math

def solve():
    """
    Finds the number of integers x for which x^3 - 16x^2 - 72x + 1056 is a perfect square.
    """
    count = 0
    
    # We determined that solutions can only exist for x in the range [-7, 62].
    for x in range(-7, 63):
        # Calculate the value of the polynomial
        val = x**3 - 16*x**2 - 72*x + 1056
        
        # Check if the value is a perfect square
        if val >= 0:
            sqrt_val = math.isqrt(val)
            if sqrt_val * sqrt_val == val:
                count += 1
                print(f"Solution found for x = {x}:")
                # As requested, output the equation for each solution.
                # Note that for negative x, parentheses are added for clarity.
                if x < 0:
                    print(f"({x})^3 - 16*({x})^2 - 72*({x}) + 1056 = {val} = {sqrt_val}^2")
                else:
                    print(f"{x}^3 - 16*{x}^2 - 72*{x} + 1056 = {val} = {sqrt_val}^2")

    print(f"\nTotal number of integers found: {count}")
    
solve()