import math

def solve():
    """
    This script finds the smallest possible denominator of the hypotenuse of a right
    triangle with area 263, all of whose sides are rational.

    The area A of such a triangle can be expressed in terms of integer generators m and n
    and a rational scaling factor k, as A = k^2 * m * n * (m^2 - n^2).
    Given A = 263, we can set k = 1/d where d is an integer.
    This leads to the condition: m * n * (m^2 - n^2) = 263 * d^2.
    
    The script searches for the smallest integers m, n (with m > n and gcd(m,n)=1)
    that satisfy this condition.

    The hypotenuse c is then given by c = (m^2 + n^2) / d. We want to find the
    smallest possible value for its denominator after simplification.
    """
    
    # We will search for m up to a reasonable limit.
    # The first solution found is expected to yield the smallest denominator.
    search_limit_m = 4000 
    
    for m in range(2, search_limit_m):
        for n in range(1, m):
            # We are looking for pairs (m, n) that are coprime.
            if math.gcd(m, n) == 1:
                
                # Calculate the product P = m * n * (m^2 - n^2)
                # m*m can be large, so ensure we use integers that can handle it.
                m_sq = m * m
                n_sq = n * n
                product_P = m * n * (m_sq - n_sq)

                # The product must be positive.
                if product_P <= 0:
                    continue

                # Check if the product P is divisible by 263.
                if product_P % 263 == 0:
                    Q = product_P // 263
                    
                    # Check if Q is a perfect square.
                    d = math.isqrt(Q)
                    if d * d == Q:
                        # A solution is found.
                        c_numerator = m_sq + n_sq
                        
                        # Simplify the fraction for the hypotenuse c = c_numerator / d
                        common_divisor = math.gcd(c_numerator, d)
                        denominator = d // common_divisor
                        
                        print(f"Found solution with m = {m}, n = {n}.")
                        print(f"The equation for the hypotenuse c is:")
                        print(f"c = (m^2 + n^2) / d")
                        print(f"c = ({m * m} + {n * n}) / {d}")
                        print(f"c = {c_numerator} / {d}")
                        print(f"The simplified denominator is {d} / gcd({c_numerator}, {d}) = {denominator}.")
                        print(f"The smallest possible denominator is {denominator}.")
                        return
                        
    print(f"No solution found within the search limit of m < {search_limit_m}.")

solve()