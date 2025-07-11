def solve():
    """
    This script calculates the value of F(n) for two prime numbers.
    The problem simplifies to calculating 27/4 mod p for each prime p.
    """
    
    primes = [80039, 80077]
    results = []

    # The simplified value to compute is (27 * 4^-1) mod p
    numerator = 27
    denominator = 4
    
    print("The problem simplifies to calculating F(n) = (27 * 4^-1) mod p for each prime p.")
    print("-" * 30)

    for p in primes:
        # Calculate the modular multiplicative inverse of the denominator
        # using Fermat's Little Theorem: a^(p-2) mod p
        # This is equivalent to pow(a, -1, p) in Python 3.8+
        mod_inverse = pow(denominator, p - 2, p)
        
        result = (numerator * mod_inverse) % p
        results.append(str(result))
        
        print(f"For p = {p}:")
        print(f"The equation is ({numerator} * {denominator}^-1) mod {p}")
        print(f"The modular inverse of {denominator} mod {p} is {mod_inverse}")
        print(f"The result is ({numerator} * {mod_inverse}) mod {p} = {result}")
        print("-" * 30)
        
    print("The final answers for p=80039 and p=80077 separated by a comma are:")
    print(",".join(results))

solve()
<<<60036,60078>>>