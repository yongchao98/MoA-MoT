import math

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    while b:
        a, b = b, a % b
    return abs(a)

def is_perfect_square(n):
    """Checks if a number is a perfect square and returns (bool, int_sqrt)."""
    if n < 0:
        return False, -1
    if n == 0:
        return True, 0
    s = int(math.sqrt(n))
    return s * s == n, s

def find_smallest_denominator(search_limit_m=4000):
    """
    Searches for the smallest possible denominator of the hypotenuse.
    
    The problem reduces to finding integers m, n such that:
    m*n*(m*m - n*n) = 263 * q*q for some integer q.
    The hypotenuse is then c = (m*m + n*n) / q.
    The denominator of c is q / gcd(q, m*m + n*n).
    We search for the (m, n) pair that minimizes this denominator.
    """
    min_denom = float('inf')
    best_solution = None

    # Search for m and n that define a primitive Pythagorean triple
    for m in range(2, search_limit_m):
        for n in range(1, m):
            # Conditions for primitive triple generators
            if (m % 2) == (n % 2): # m and n must have opposite parity
                continue
            if gcd(m, n) != 1: # m and n must be coprime
                continue

            # This is the area of the integer Pythagorean triple (m^2-n^2, 2mn, m^2+n^2)
            # which is Area_int = m * n * (m*m - n*n)
            # In our case, Area_int = 263 * q^2
            area_int = m * n * (m*m - n*n)
            
            if area_int % 263 == 0:
                quotient = area_int // 263
                is_sq, q = is_perfect_square(quotient)

                if is_sq:
                    hypotenuse_numerator = m*m + n*n
                    
                    common_divisor = gcd(q, hypotenuse_numerator)
                    denominator = q // common_divisor

                    if denominator < min_denom:
                        min_denom = denominator
                        best_solution = (m, n, q, denominator)

    if best_solution:
        m, n, q, denom = best_solution
        print(f"Smallest denominator found up to m={search_limit_m}: {min_denom}")
        print(f"Achieved with m = {m}, n = {n}")
        a_num = m*m - n*n
        b_num = 2*m*n
        c_num = m*m + n*n
        
        # We need to simplify the fractions to their lowest terms
        # Sides are (a_num/q), (b_num/q). Hypotenuse is (c_num/q).
        
        common_a = gcd(a_num, q)
        a_side = f"{a_num//common_a}/{q//common_a}"

        common_b = gcd(b_num, q)
        b_side = f"{b_num//common_b}/{q//common_b}"
        
        common_c = gcd(c_num, q)
        c_side = f"{c_num//common_c}/{q//common_c}"

        print("\nCorresponding rational triangle:")
        print(f"Side a = {a_side}")
        print(f"Side b = {b_side}")
        print(f"Hypotenuse c = {c_side}")
        
    else:
        print(f"No solution found up to m={search_limit_m}. You might need to increase the search limit.")

if __name__ == '__main__':
    # The solution involves large m and n, so a higher search limit is needed
    # than typical for simple search problems. Let's try up to a few thousand.
    find_smallest_denominator(search_limit_m=4000)
