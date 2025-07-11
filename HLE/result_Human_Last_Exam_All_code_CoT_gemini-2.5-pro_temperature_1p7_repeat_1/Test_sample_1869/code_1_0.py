def find_special_prime():
    """
    Finds the smallest prime p from a given list such that
    Z[p-th root of 6] is not the ring of integers of Q(p-th root of 6).

    This is equivalent to finding the smallest prime p in the list for which
    6^(p-1) is congruent to 1 modulo p^2.
    """
    # The list of primes from the answer choices
    primes = [17, 383, 1093, 66161, 534851]
    base = 6

    print("We are looking for the smallest prime p > 3 in the choices such that Z[p-th root of 6] is not the ring of integers of Q(p-th root of 6).")
    print("This condition is equivalent to finding the smallest prime p such that 6^(p-1) is congruent to 1 modulo p^2.")
    print("-" * 70)

    for p in primes:
        p_squared = p * p
        exponent = p - 1
        
        print(f"Checking prime p = {p}...")
        
        # Using Python's built-in pow for efficient modular exponentiation
        result = pow(base, exponent, p_squared)
        
        print(f"We compute: {base}^({p} - 1) mod {p}^2")
        print(f"Result: {base}^{exponent} mod {p_squared} = {result}")
        
        if result == 1:
            print(f"\nThe condition 6^(p-1) = 1 (mod p^2) is met for p = {p}.")
            print(f"This is the smallest prime in the list that satisfies the condition.")
            print("\nThe final equation with each number explicitly stated is:")
            print(f"  Base: {base}")
            print(f"  Prime p: {p}")
            print(f"  Exponent (p - 1): {exponent}")
            print(f"  Modulus (p^2): {p_squared}")
            print(f"The equation {base}^{exponent} \u2261 1 (mod {p_squared}) holds.")
            return
        else:
            print(f"The result is not 1. The condition is not met for p = {p}.")
            print("-" * 70)

find_special_prime()