import math

def solve():
    """
    This function counts the number of distinct parallelograms according to the given restrictions.
    """
    count = 0
    limit = 100
    
    # We assume a < b to avoid counting (a,b) and (b,a) as different side pairs.
    for a in range(1, limit):
        # The condition a + b < 100 implies b < 100 - a.
        for b in range(a + 1, limit - a):
            
            # Condition 2: a and b are coprime integers.
            if math.gcd(a, b) != 1:
                continue

            # Condition 1 (a!=b) and Condition 2 (2a < a+b) are guaranteed by the loop ranges.
            
            # From the parallelogram law, d1^2 + d2^2 = 2(a^2 + b^2).
            # This is our target number N for the sum of two squares.
            N = 2 * (a**2 + b**2)
            
            # Find integer diagonals d1, d2. We search for d1 and calculate the corresponding d2.
            # To avoid duplicate diagonal pairs {d1,d2}, we enforce d1 < d2.
            # If d1^2 + d2^2 = N, then d1^2 must be less than N/2.
            d1_limit = int(math.sqrt(N / 2))
            for d1 in range(1, d1_limit + 1):
                d2_squared = N - d1**2
                d2_sqrt_int = int(math.sqrt(d2_squared))

                if d2_sqrt_int * d2_sqrt_int == d2_squared:
                    d2 = d2_sqrt_int
                    # At this point, we have integer diagonals d1 and d2 where d1 < d2.
                    # This satisfies Condition 4 and the 'not a rectangle' part of Condition 1.
                    
                    # Now, check for integer area (Condition 3).
                    # The parallelogram's area is integer if the triangle (a, b, d1) is Heronian
                    # and its area is a multiple of 1/2.
                    # 16 * Area_triangle^2 = ((a+b)^2 - d1^2) * (d1^2 - (b-a)^2)
                    
                    P = ((a + b)**2 - d1**2) * (d1**2 - (b - a)**2)
                    
                    # P must be positive for a non-degenerate triangle.
                    if P > 0:
                        S_int = int(math.sqrt(P))
                        # P must be a perfect square for the area to be rational.
                        if S_int * S_int == P:
                            # Area_parallelogram = sqrt(P) / 2.
                            # For this to be an integer, S must be an even integer.
                            if S_int % 2 == 0:
                                # All conditions are met. This is a valid parallelogram.
                                count += 1
                                
    print(count)

solve()