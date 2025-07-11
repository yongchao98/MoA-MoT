import math

def solve():
    """
    This function calculates the number of distinct parallelograms satisfying the given conditions.
    
    The conditions are:
    1. It's not a rhombus (a != b) and not a rectangle (d1 != d2).
    2. Sides a, b are coprime integers with a < b and a + b < 100.
    3. The area is an integer.
    4. Both diagonals d1, d2 are integers.
    """
    
    count = 0
    
    # Iterate through all possible side lengths a and b
    # From 2a < a + b < 100, we get a < b and a + b < 100.
    # The smallest a can be is 1.
    # The largest a can be is 49 (e.g., a=49, b=50 gives a+b=99).
    for a in range(1, 50):
        # b must be greater than a and less than 100 - a
        for b in range(a + 1, 100 - a):
            
            # Condition 2: a and b must be coprime.
            if math.gcd(a, b) != 1:
                continue

            # Parallelogram law: d1^2 + d2^2 = 2*(a^2 + b^2)
            s = 2 * (a**2 + b**2)

            # Find integer diagonals d1, d2.
            # To ensure d1 < d2, we loop d1 up to sqrt(s/2).
            # We exclude d1 = d2, which would be a rectangle.
            d1_limit = math.isqrt(s // 2)
            for d1 in range(b - a + 1, d1_limit + 1):
                
                d2_squared = s - d1**2
                
                # We need d1 < d2, so d1^2 < d2_squared
                if d1**2 >= d2_squared:
                    continue

                d2 = math.isqrt(d2_squared)

                # Check if d2 is an integer
                if d2**2 == d2_squared:
                    
                    # We have integer diagonals (d1, d2)
                    # The loops already ensure a!=b (rhombus) and d1!=d2 (rectangle).

                    # Check triangle inequalities for non-degenerate parallelograms.
                    # This means |a-b| < d1 < a+b and |a-b| < d2 < a+b.
                    # Since we iterate d1 > b-a, we just need to check d2 < a+b.
                    if not (d2 < a + b):
                        continue

                    # Check for integer area.
                    # Area^2 = ((a+b)^2 - d1^2) * (d1^2 - (b-a)^2) / 4.
                    # So, the numerator must be a perfect square and its sqrt must be even.
                    val = ((a + b)**2 - d1**2) * (d1**2 - (b - a)**2)
                    
                    if val <= 0:
                        continue # Degenerate case

                    sqrt_val = math.isqrt(val)
                    if sqrt_val**2 == val:
                        # For area to be integer, sqrt_val must be divisible by 2.
                        if sqrt_val % 2 == 0:
                            count += 1
    
    print(f"The number of distinct parallelograms is:")
    print(count)

solve()