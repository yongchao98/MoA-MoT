def find_special_prime():
    """
    This function finds the smallest prime p from a given list of choices
    such that Z[6^(1/p)] is not the ring of integers of Q(6^(1/p)).
    This is equivalent to finding the smallest p such that 6^(p-1) = 1 (mod p^2).
    """

    # Answer choices for p, in ascending order.
    primes = [17, 383, 1093, 66161, 534851]
    base = 6
    found_prime = None

    print("Searching for the smallest prime p in the list such that 6^(p-1) is congruent to 1 mod p^2.\n")

    for p in primes:
        modulus = p * p
        exponent = p - 1

        # Use Python's pow(base, exp, mod) for efficient modular exponentiation.
        result = pow(base, exponent, modulus)

        print(f"Checking prime p = {p}...")
        
        if result == 1:
            found_prime = p
            print(f"Condition met for p = {p}.")
            print("The final equation is:")
            print(f"{base}^({p} - 1) ≡ 1 (mod {p}^2)")
            print(f"Substituting the values:")
            print(f"{base}^{exponent} ≡ 1 (mod {modulus})")
            print(f"The calculation pow({base}, {exponent}, {modulus}) results in {result}, confirming the condition.")
            break # Stop after finding the smallest prime in the list.
        else:
            print(f"For p = {p}, 6^({p-1}) mod {p*p} = {result}. Condition not met.")
            print("-" * 30)

    if found_prime is None:
        print("\nNo prime in the list satisfies the condition.")

find_special_prime()