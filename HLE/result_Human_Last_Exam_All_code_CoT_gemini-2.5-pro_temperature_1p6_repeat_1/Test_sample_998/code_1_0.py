import collections

def get_prime_factorization(n):
    """
    Returns the prime factorization of n as a collections.Counter.
    """
    factors = collections.Counter()
    d = 2
    temp_n = n
    while d * d <= temp_n:
        while temp_n % d == 0:
            factors[d] += 1
            temp_n //= d
        d += 1
    if temp_n > 1:
        factors[temp_n] += 1
    return factors

def solve():
    """
    Solves the problem as described in the plan.
    """
    # Step 2: Define exponents from N = 2^8 * 29^59 * 59^79 * 79^29
    a = 8  # for prime 2
    b = 59 # for prime 29 (1 mod 4)
    c1 = 79 # for prime 59 (3 mod 4)
    c2 = 29 # for prime 79 (3 mod 4)

    # Step 3: Evaluate the components of S(N) = term_a * P1 * P2
    term_a = a + 1
    
    # P1 is the sum related to the prime congruent to 1 mod 4
    p1 = (b + 1) * (b + 2) // 2

    # P2 is the sum related to primes congruent to 3 mod 4
    # S_e(c) = sum_{k=0, even}^c (k+1) = (floor(c/2)+1)^2
    # S_o(c) = sum_{k=0, odd}^c (k+1) = (floor((c-1)/2)+1)*(floor((c-1)/2)+2)
    # For c1 = 79, m1 = 39
    m1 = c1 // 2
    se_c1 = (m1 + 1)**2
    so_c1 = (m1 + 1) * (m1 + 2)
    
    # For c2 = 29, m2 = 14
    m2 = c2 // 2
    se_c2 = (m2 + 1)**2
    so_c2 = (m2 + 1) * (m2 + 2)

    # c1 and c2 have the same parity (both odd), so C1 and C2 must have same parity
    p2 = se_c1 * se_c2 + so_c1 * so_c2
    
    # Step 4: Find the prime factorization of S(N)
    factors_a = get_prime_factorization(term_a)
    factors_p1 = get_prime_factorization(p1)
    factors_p2 = get_prime_factorization(p2)
    
    total_factors = factors_a + factors_p1 + factors_p2
    
    # Calculate the number of divisors from the combined factors
    num_divisors = 1
    
    # Build the expression string for the final output
    expr_parts = []
    
    # Sort keys for consistent output order
    sorted_primes = sorted(total_factors.keys())
    
    for prime in sorted_primes:
        exponent = total_factors[prime]
        num_divisors *= (exponent + 1)
        expr_parts.append(f"({exponent}+1)")
        
    final_equation = " * ".join(expr_parts)
    print(f"The number of divisors is calculated from the exponents of the prime factors of S(N).")
    print(f"S(N) = {term_a} * {p1} * {p2}")
    prime_factorization_str = " * ".join([f"{p}^{e}" for p, e in sorted(total_factors.items())])
    print(f"The prime factorization of S(N) is {prime_factorization_str}")
    print(f"Number of divisors = {final_equation} = {num_divisors}")
    
    return num_divisors

# Run the solver and capture the final numerical answer.
final_answer = solve()
print(f"<<<{final_answer}>>>")