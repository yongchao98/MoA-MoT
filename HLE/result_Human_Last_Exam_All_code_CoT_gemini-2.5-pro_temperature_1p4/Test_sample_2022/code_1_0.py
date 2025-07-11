def power(a, b, m):
    """Computes (a^b) % m using modular exponentiation."""
    res = 1
    a %= m
    while b > 0:
        if b % 2 == 1:
            res = (res * a) % m
        a = (a * a) % m
        b //= 2
    return res

def inv(n, m):
    """Computes modular inverse of n modulo m using Fermat's Little Theorem."""
    return power(n, m - 2, m)

def solve():
    """
    Solves the problem for the two given primes.
    """
    primes = [80039, 80077]
    results = []

    # Case p = 80039
    p1 = primes[0]
    # We determined F(N) = S(3) = 7/8 mod p1
    # We print the calculation for S(3)
    s0, s1, s2 = 1, 1, 1
    s3_num = 8*s2 - s0
    s3_den = 8
    
    val1 = (s3_num * inv(s3_den, p1)) % p1
    results.append(val1)

    print("For p = 80039:")
    print(f"The argument N reduces to 3 modulo (p-1).")
    print("F(N) is congruent to S(3) mod p.")
    print("S(0) = 1")
    print("S(1) = 1")
    print("S(2) = 1")
    print("The recurrence is S(n) = S(n-1) - (1/8) * S(n-3), so 8*S(n) = 8*S(n-1) - S(n-3).")
    print(f"8*S(3) = 8*S(2) - S(0) = 8*({s2}) - {s0} = {s3_num}")
    print(f"S(3) = {s3_num}/{s3_den}")
    print(f"F(N) = {s3_num} * ({s3_den}^-1 mod {p1}) mod {p1} = {s3_num} * {inv(s3_den, p1)} mod {p1} = {val1}")
    print("-" * 20)
    
    # Case p = 80077
    p2 = primes[1]
    # We determined F(N) = S(7) = 27/64 mod p2
    # We print the calculation for S(7)
    s = [1, 1, 1]
    s_num_den = [(8,8), (8,8), (8,8)] # Store S(n) as fractions a/b -> (8a, 8b) for common denominator
    
    # S(3) = 7/8
    s.append(7/8)
    s_num_den.append((7, 8))

    # S(4) = S(3) - S(1)/8 = 7/8 - 1/8 = 6/8
    s.append(6/8)
    s_num_den.append((6, 8))
    
    # S(5) = S(4) - S(2)/8 = 6/8 - 1/8 = 5/8
    s.append(5/8)
    s_num_den.append((5, 8))

    # S(6) = S(5) - S(3)/8 = 5/8 - (7/8)/8 = 40/64 - 7/64 = 33/64
    s.append(33/64)
    s_num_den.append((33, 64))

    # S(7) = S(6) - S(4)/8 = 33/64 - (6/8)/8 = 33/64 - 6/64 = 27/64
    s.append(27/64)
    s_num_den.append((27, 64))
    
    s7_num = 27
    s7_den = 64
    val2 = (s7_num * inv(s7_den, p2)) % p2
    results.append(val2)
    
    print("For p = 80077:")
    print(f"The argument N reduces to 7 modulo (p+1).")
    print("We proceed with the hypothesis that F(N) is congruent to S(7) mod p.")
    print("S(4) = S(3) - S(1)/8 = 7/8 - 1/8 = 6/8")
    print("S(5) = S(4) - S(2)/8 = 6/8 - 1/8 = 5/8")
    print("S(6) = S(5) - S(3)/8 = 5/8 - (7/8)/8 = 33/64")
    print("S(7) = S(6) - S(4)/8 = 33/64 - (6/8)/8 = 27/64")
    print(f"F(N) = {s7_num} * ({s7_den}^-1 mod {p2}) mod {p2} = {s7_num} * {inv(s7_den, p2)} mod {p2} = {val2}")
    print("-" * 20)

    print("Final answer:")
    print(','.join(map(str, results)))

solve()