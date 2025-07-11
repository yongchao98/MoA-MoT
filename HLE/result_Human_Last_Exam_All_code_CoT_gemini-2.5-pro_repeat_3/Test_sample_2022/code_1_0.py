def solve():
    """
    This function calculates the value of F(n) for two prime numbers.
    The formula for F(n) is derived from a recurrence relation and properties of sequences over finite fields.
    The final value is F(3) = 1 - 8^(-1) mod p.
    """
    
    # First prime number
    p1 = 80039
    
    # Calculate modular inverse of 8 mod p1
    inv8_1 = pow(8, -1, p1)
    
    # Calculate F(3) for p1
    # F(3) = 1 - 8^(-1) mod p1
    ans1 = (1 - inv8_1 + p1) % p1
    
    # Second prime number
    p2 = 80077
    
    # Calculate modular inverse of 8 mod p2
    inv8_2 = pow(8, -1, p2)
    
    # Calculate F(3) for p2
    # F(3) = 1 - 8^(-1) mod p2
    ans2 = (1 - inv8_2 + p2) % p2
    
    # The problem asks for the calculated values of F(n) for the given primes.
    # Our derivation shows that F(n) simplifies to F(3) in both cases.
    # Let's show the value of n modulo p for context, although it's not strictly needed for the calculation.
    n_mod_p1 = (p1**5 + 2*p1**4 - 19*p1**3 - 3*p1**2 + 16*p1 + 6) % p1
    n_mod_p2 = (p2**5 + 2*p2**4 - 19*p2**3 - 3*p2**2 + 16*p2 + 6) % p2

    # The recurrence relation is F(n) = F(n-1) - (8^(-1) mod p) * F(n-3)
    # The argument n reduces to 3.
    # We need to compute F(3) = F(2) - (8^(-1) mod p) * F(0)
    # With F(0) = 1, F(1) = 1, F(2) = 1
    # F(3) = 1 - 8^(-1) mod p
    
    print(f"{ans1},{ans2}")

solve()