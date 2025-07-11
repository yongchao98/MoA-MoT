def solve():
    """
    Calculates the number of divisors of S, based on the derived formula.
    S = 9 * (sum over b) * (sum over c,d)
    """

    # Prime exponents of N = 2^e2 * 29^e29 * 59^e59 * 79^e79
    e2 = 8
    e29 = 59
    e59 = 79
    e79 = 29

    # Constant term from tau(2^8)
    const_term = e2 + 1
    # print(f"Constant term: {const_term}")

    # First sum: sum_{b=0}^{59} (60-b)
    # This is an arithmetic series sum 1 + 2 + ... + 60
    sum_b_val = (e29 + 1) * (e29 + 2) // 2
    # print(f"Sum over b: {sum_b_val}")

    # Second sum: sum over c and d
    # Sum over c (even)
    sum_c_even = sum(e59 + 1 - c for c in range(0, e59 + 1, 2))
    # Sum over c (odd)
    sum_c_odd = sum(e59 + 1 - c for c in range(1, e59 + 1, 2))
    # Sum over d (even)
    sum_d_even = sum(e79 + 1 - d for d in range(0, e79 + 1, 2))
    # Sum over d (odd)
    sum_d_odd = sum(e79 + 1 - d for d in range(1, e79 + 1, 2))
    
    sum_cd_val = sum_c_even * sum_d_even + sum_c_odd * sum_d_odd
    # print(f"Sum over c,d: {sum_cd_val}")
    
    # We found the prime factorization of S to be:
    # S = 2^7 * 3^4 * 5^3 * 61^1 * 157^1
    # So the exponents are e_i = {7, 4, 3, 1, 1}
    
    exponents = [7, 4, 3, 1, 1]
    
    num_divisors = 1
    for exp in exponents:
        num_divisors *= (exp + 1)
        
    factors_str = " * ".join([f"({e}+1)" for e in exponents])
    print(f"The number of divisors is the product of (exponent + 1) for each prime factor of S.")
    print(f"S has prime exponents: {exponents}")
    print(f"Number of divisors = {factors_str} = {num_divisors}")

solve()
<<<640>>>