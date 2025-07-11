def find_prime():
    """
    Finds the smallest prime p in a given list such that 6^(p-1) is congruent to 1 modulo p^2.
    This is the condition for Z[p-th root of 6] not being the ring of integers of Q(p-th root of 6).
    """
    # Answer choices provided in the problem
    primes = [17, 383, 1093, 66161, 534851]
    
    print("The condition for Z[p-th root of 6] not to be the ring of integers is 6^(p-1) ≡ 1 (mod p^2).")
    print(f"Checking the primes from the list: {primes}\n")

    found_prime = None
    
    for p in primes:
        # Calculate 6^(p-1) mod p^2
        # pow(base, exp, mod) is efficient for modular exponentiation
        modulus = p * p
        result = pow(6, p - 1, modulus)
        
        print(f"Checking for p = {p}:")
        print(f"  6^({p-1}) mod {modulus} = {result}")
        
        if result == 1:
            found_prime = p
            print(f"\nCondition met for p = {p}.")
            break
        else:
            print(f"  Condition not met for p = {p}.\n")

    if found_prime:
        p = found_prime
        base = 6
        exponent = p - 1
        modulus = p * p
        result = pow(base, exponent, modulus)
        
        print("\nThe smallest prime in the list that satisfies the condition is", p)
        print("The final equation is:")
        print(f"{base}^{exponent} mod {modulus} = {result}")

find_prime()