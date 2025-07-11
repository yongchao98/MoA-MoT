import sys

def find_special_prime():
    """
    Finds the smallest prime p > 3 from a list of choices such that
    Z[p√6] is not the ring of integers of Q(p√6).

    This is equivalent to finding the smallest prime p in the list
    that satisfies the congruence 6^(p-1) ≡ 1 (mod p^2).
    """
    # The given answer choices
    primes = [17, 383, 1093, 66161, 534851]

    found_prime = None

    # Iterate through the primes in ascending order to find the smallest one
    for p in primes:
        p_squared = p * p
        # Check the condition 6^(p-1) ≡ 1 (mod p^2)
        if pow(6, p - 1, p_squared) == 1:
            found_prime = p
            break
    
    if found_prime is None:
        print("None of the primes in the list satisfy the condition.")
        return

    # Unpack the variables for the final equation
    p = found_prime
    base = 6
    exponent = p - 1
    modulus = p * p
    result = 1 # We know this from the check that found the prime

    print(f"The condition for Z[p√6] to not be the ring of integers is 6^(p-1) ≡ 1 (mod p^2).")
    print(f"Checking the primes {primes}, we find the smallest prime satisfying this is p = {p}.")
    print("\nThe final equation is:")
    print(f"6^{exponent} ≡ 1 (mod {modulus})")
    print("\nThe numbers in this equation are:")
    print(f"base = {base}")
    print(f"exponent = {exponent}")
    print(f"modulus = {modulus}")
    print(f"result = {result}")

find_special_prime()
