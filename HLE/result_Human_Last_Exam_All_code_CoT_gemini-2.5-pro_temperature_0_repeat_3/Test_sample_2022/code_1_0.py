def solve():
    """
    Calculates the value of F(n) for two prime numbers.
    The problem simplifies to calculating S(6) mod p.
    S(6) = 33/64.
    We calculate (33 * 64^-1) mod p for each prime.
    """
    primes = [80039, 80077]
    results = []
    
    # The value of S(6) is 33/64
    a = 33
    b = 64

    print(f"The expression for n is p^5+2p^4-19p^3-3p^2+16p+6.")
    print(f"Modulo p, n is congruent to 6.")
    print(f"We calculate S(6) = sum_{k=0 to 3} C(6-2k, k) * (-1/8)^k")
    print(f"S(6) = C(6,0)*(-1/8)^0 + C(4,1)*(-1/8)^1 + C(2,2)*(-1/8)^2")
    print(f"S(6) = 1*1 + 4*(-1/8) + 1*(1/64) = 1 - 1/2 + 1/64 = 33/64")
    print("-" * 20)

    for p in primes:
        print(f"Calculating for p = {p}:")
        
        # Calculate modular inverse of b modulo p
        # We need to solve (b * x) % p = 1
        try:
            b_inv = pow(b, -1, p)
        except ValueError:
            print(f"Error: {b} is not invertible modulo {p}")
            continue

        # Calculate F(n) = (a * b_inv) mod p
        result = (a * b_inv) % p
        results.append(result)
        
        print(f"F(n) mod {p} = S(6) mod {p} = (33 * 64^-1) mod {p}")
        print(f"The value of 64^-1 mod {p} is: {b_inv}")
        print(f"The equation is: ({a} * {b_inv}) mod {p}")
        print(f"The result is: {result}")
        print("-" * 20)

    print(f"Final answers separated by a comma:")
    print(','.join(map(str, results)))

solve()