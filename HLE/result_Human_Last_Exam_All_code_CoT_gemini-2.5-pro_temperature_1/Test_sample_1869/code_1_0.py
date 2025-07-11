def find_special_prime():
    """
    Finds the smallest prime p from a list such that the ring Z[p-th root of 6]
    is not the ring of integers of Q(p-th root of 6).

    This is equivalent to finding the smallest prime p in the list
    that satisfies the condition: 6^(p-1) = 1 (mod p^2).
    """
    primes = [17, 383, 1093, 66161, 534851]
    
    print("We are looking for the smallest prime p such that 6^(p-1) = 1 (mod p^2).")
    print("We will test the primes from the answer choices in order.")
    print("-" * 40)
    
    found_prime = None
    
    for p in primes:
        p_squared = p * p
        
        # Use Python's built-in pow() for efficient modular exponentiation.
        # pow(base, exponent, modulus)
        result = pow(6, p - 1, p_squared)
        
        print(f"Checking for p = {p}:")
        print(f"Is 6^({p-1}) = 1 (mod {p_squared})?")
        
        if result == 1:
            print(f"Result: {result}. The condition is satisfied.")
            found_prime = p
            break
        else:
            print(f"Result: {result}. The condition is not satisfied.")
            print("-" * 40)

    if found_prime:
        p = found_prime
        p_squared = p * p
        p_minus_1 = p - 1
        print("\n" + "=" * 40)
        print(f"The smallest prime in the list that satisfies the condition is {p}.")
        print("The final equation is:")
        print(f"6^({p} - 1) = 1 (mod {p}^2)")
        print("Substituting the values:")
        print(f"6^{p_minus_1} = 1 (mod {p_squared})")

find_special_prime()