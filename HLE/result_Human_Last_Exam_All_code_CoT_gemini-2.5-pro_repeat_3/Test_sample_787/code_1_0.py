def solve():
    """
    This function calculates the limit of the sequence g_n.
    The limit L is given by the greatest common divisor of P(p) for all sufficiently large primes p.
    This implies L is the product of q^(v_q(L)) for all primes q.
    
    The prime factors of L can only be 2, 3, and 5.
    
    The exponents (p-adic valuations) are calculated as:
    v_q(L) = min_{a in Z, q does not divide a} v_q(P(a))
    
    - For q=2, the exponent is 10.
    - For q=3, the exponent is 2.
    - For q=5, the exponent is 1.
    - For all other primes q, the exponent is 0.
    
    So, L = 2^10 * 3^2 * 5^1.
    """
    
    p1 = 2
    e1 = 10
    
    p2 = 3
    e2 = 2
    
    p3 = 5
    e3 = 1
    
    val1 = p1 ** e1
    val2 = p2 ** e2
    val3 = p3 ** e3
    
    limit = val1 * val2 * val3
    
    print(f"{p1}^{e1} * {p2}^{e2} * {p3}^{e3} = {limit}")

solve()
<<<46080>>>