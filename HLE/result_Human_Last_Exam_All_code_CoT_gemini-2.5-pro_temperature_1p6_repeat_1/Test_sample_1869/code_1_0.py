def find_special_prime():
    """
    This function finds the smallest prime p from a given list such that
    the ring Z[6^(1/p)] is not the ring of integers of Q(6^(1/p)).

    This condition is equivalent to checking if 6^(p-1) is congruent to 1 modulo p^2.
    """
    # Primes from the answer choices in increasing order
    primes = [17, 383, 1093, 66161, 534851]

    print("We are looking for the smallest prime p such that 6^(p-1) = 1 (mod p^2).\n")

    for p in primes:
        print(f"--- Checking prime p = {p} ---")
        p_squared = p * p
        exponent = p - 1
        base = 6

        # Calculate 6^(p-1) mod p^2 using modular exponentiation
        print(f"Calculating: {base}^({exponent}) mod {p_squared}")
        result = pow(base, exponent, p_squared)
        print(f"Result: {result}")

        if result == 1:
            print(f"\nThe condition 6^({p-1}) \u2261 1 (mod {p*p}) is met for p = {p}.")
            print("This is the smallest prime in the list that satisfies the condition.")
            
            # As requested, printing the numbers in the final equation
            print("\nThe numbers in the final equation are:")
            print(f"Base: {base}")
            print(f"Exponent: {exponent}")
            print(f"Result of congruence: 1")
            print(f"Modulus: {p_squared}")
            
            return p
        else:
            print(f"The result is not 1, so this is not the answer. Moving to the next prime.\n")

# Run the function to find the answer
find_special_prime()