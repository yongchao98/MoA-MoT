def find_special_prime():
    """
    This function finds the smallest prime p from a list of choices such that
    Z[p-th root of 6] is not the ring of integers of Q(p-th root of 6).
    The condition for this is 6^(p-1) is congruent to 1 modulo p^2.
    """
    
    primes = [17, 383, 1093, 66161, 534851]
    
    print("We are looking for the smallest prime p such that 6^(p-1) == 1 (mod p^2).")
    print("Let's check the primes from the answer choices:")
    print("-" * 40)

    for p in primes:
        # Calculate p^2
        p_squared = p * p
        
        # Calculate 6^(p-1) mod p^2 using Python's efficient pow() function
        result = pow(6, p - 1, p_squared)
        
        # Print the equation for the current prime p
        print(f"For p = {p}, the calculation is:")
        print(f"6^({p}-1) mod {p}^2 = {result}")
        
        if result == 1:
            print(f"\nThe condition is met for p = {p}.")
            print(f"Therefore, {p} is the smallest prime in the list for which Z[p-th root of 6] is not the ring of integers.")
            break
        else:
            print("The condition is not met.")
            print("-" * 40)

find_special_prime()