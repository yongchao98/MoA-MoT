def find_special_prime():
    """
    Finds the smallest prime p from a given list of choices such that
    the ring Z[p-th root of 6] is not the ring of integers of Q(p-th root of 6).

    This is equivalent to finding the smallest prime p in the list
    that satisfies the condition: 6^(p-1) is congruent to 1 modulo p^2.
    """
    
    # Answer choices given as a list of primes
    primes = [17, 383, 1093, 66161, 534851]
    a = 6

    # Iterate through the primes in the order they are given (which is ascending)
    for p in primes:
        # Calculate p^2 for the modulus
        p_squared = p * p
        
        # Calculate a^(p-1) mod p^2
        # Use pow(base, exp, mod) for efficient computation of modular exponentiation
        result = pow(a, p - 1, p_squared)
        
        # Check if the condition is met
        if result == 1:
            print(f"Found the smallest prime p = {p} in the list that satisfies the condition.")
            print(f"The condition is that a^(p-1) is congruent to 1 modulo p^2.")
            print("\nLet's verify the equation for p = {p}:")
            print(f"a^(p-1) \u2261 1 (mod p^2)")
            print(f"{a}^({p}-1) \u2261 1 (mod {p}^2)")
            print(f"{a}^{p-1} \u2261 1 (mod {p_squared})")
            
            print("\nNumbers in the final equation:")
            print(f"Base a = {a}")
            print(f"Prime p = {p}")
            print(f"Exponent (p-1) = {p-1}")
            print(f"Modulus (p^2) = {p_squared}")
            print(f"Calculation: pow({a}, {p-1}, {p_squared}) returns {result}.")
            print(f"\nSince the result is 1, p = {p} is the answer.")
            
            # Stop after finding the first (and smallest) prime that satisfies the condition
            return

# Execute the function
find_special_prime()