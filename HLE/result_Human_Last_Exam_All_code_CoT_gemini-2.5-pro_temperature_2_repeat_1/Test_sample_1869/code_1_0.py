def find_irregular_prime():
    """
    This function searches for the smallest prime p from a given list
    such that 6^(p-1) is congruent to 1 modulo p^2. This condition
    determines if Z[p-th root of 6] is the ring of integers of
    Q(p-th root of 6).
    """
    
    primes = [17, 383, 1093, 66161, 534851]
    
    print("Finding the smallest prime p such that 6^(p-1) = 1 (mod p^2)")
    print("-" * 60)

    found_prime = None
    
    for p in primes:
        base = 6
        exponent = p - 1
        modulus = p * p
        
        # Calculate 6^(p-1) mod p^2
        result = pow(base, exponent, modulus)
        
        print(f"Checking for p = {p}:")
        print(f"Equation: {base}^({p}-1) mod {p}^2")
        print(f"Calculation: {base}^{exponent} mod {modulus} = {result}")
        
        if result == 1:
            found_prime = p
            print(f"\nCondition is met for p = {p}.")
            print(f"The smallest prime in the list is {p}.")
            break
        else:
            print("Condition not met.\n")
            
    if not found_prime:
        print("No prime in the list satisfies the condition.")

if __name__ == "__main__":
    find_irregular_prime()