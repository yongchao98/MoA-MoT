def p_adic_valuation(n, p):
    """Calculates the p-adic valuation of n, i.e., the exponent of p in the prime factorization of n."""
    if n == 0:
        return float('inf') 
    count = 0
    # We can only divide by p if n is a multiple of p
    while n > 0 and n % p == 0:
        count += 1
        n //= p
    return count

def find_smallest_n():
    """
    Finds the smallest integer n >= 2 satisfying the two properties.
    The properties translate to two possible cases for n:
    Case A: n is odd, v_2(n-1)=9, and n is a multiple of 5.
    Case B: n is even, v_5(n-1)=9, and n is not a multiple of 5.
    """
    n = 2
    while True:
        n_minus_1 = n - 1

        v2_n = p_adic_valuation(n, 2)
        v5_n = p_adic_valuation(n, 5)
        
        # In Case A, v_p(n-1) depends on the other primes.
        v2_n_minus_1 = p_adic_valuation(n_minus_1, 2)
        # In Case B, v_p(n-1) depends on the other primes.
        v5_n_minus_1 = p_adic_valuation(n_minus_1, 5)
        
        # Case A: n is odd (v2_n=0), v_2(n-1) is exactly 9, and n is a multiple of 5 (v5_n>0)
        cond_A = (v2_n == 0 and v2_n_minus_1 == 9 and v5_n > 0)
        
        # Case B: n is not a multiple of 5 (v5_n=0), v_5(n-1) is exactly 9, and n is even (v2_n>0)
        cond_B = (v5_n == 0 and v5_n_minus_1 == 9 and v2_n > 0)
        
        if cond_A or cond_B:
            return n
        
        n += 1

result = find_smallest_n()
print(f"The smallest positive integer n is: {result}")