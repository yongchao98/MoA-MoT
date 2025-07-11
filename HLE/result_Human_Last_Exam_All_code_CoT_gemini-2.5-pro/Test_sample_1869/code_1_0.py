def find_special_prime():
    """
    Finds the smallest prime p from the given choices such that
    Z[p-th root of 6] is not the ring of integers of Q(p-th root of 6).

    This occurs when 6^(p-1) is congruent to 1 modulo p^2.
    """
    
    # Answer choices given in the problem
    answer_choices = {
        'A': 17,
        'B': 383,
        'C': 1093,
        'D': 66161,
        'E': 534851
    }
    
    # We test the primes in ascending order to find the smallest one.
    primes_to_test = sorted(answer_choices.values())
    
    print("Searching for the smallest prime p such that 6^(p-1) = 1 (mod p^2)...")
    print("-" * 60)

    found_prime = None
    found_letter = None

    for p in primes_to_test:
        p_squared = p * p
        # Calculate 6^(p-1) mod p^2
        result = pow(6, p - 1, p_squared)
        
        print(f"Testing prime p = {p}...")
        
        if result == 1:
            print(f"Condition met for p = {p}!")
            print("This means Z[p-th root of 6] is NOT the ring of integers.")
            
            # Find the corresponding letter for the choice
            for letter, value in answer_choices.items():
                if value == p:
                    found_letter = letter
                    break
            
            found_prime = p
            
            # Print the final equation with all its components
            print("\nThe final equation is:")
            base = 6
            exponent = p - 1
            modulus = p_squared
            
            print(f"{base}^{{{exponent}}} \u2261 {result} (mod {modulus})")
            print(f"Where:")
            print(f"  Base     = {base}")
            print(f"  Exponent = p - 1 = {p} - 1 = {exponent}")
            print(f"  Modulus  = p^2 = {p}^2 = {modulus}")
            print(f"  Result   = {result}")

            break # Stop after finding the smallest prime
        else:
            print(f"Condition not met. 6^({p}-1) mod {p_squared} = {result}, not 1.")
            print("-" * 60)

    if found_prime is None:
        print("\nNo prime in the list satisfies the condition.")
    else:
        print(f"\nThe smallest prime from the choices is {found_prime}, which is choice {found_letter}.")
        print(f"<<<{found_letter}>>>")

# Execute the function
find_special_prime()