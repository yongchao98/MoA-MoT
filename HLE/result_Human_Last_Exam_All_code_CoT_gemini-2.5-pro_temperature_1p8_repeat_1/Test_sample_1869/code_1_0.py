import sys

def solve():
    """
    Finds the smallest prime p > 3 from a list such that Z[6^(1/p)] is not the ring of integers of Q(6^(1/p)).
    """
    # This problem from algebraic number theory reduces to a condition in modular arithmetic.
    # The ring Z[6^(1/p)] is not the ring of integers if and only if the prime p
    # satisfies the congruence 6^(p-1) ≡ 1 (mod p^2).
    # Such a prime is called a Wieferich prime to base 6.

    # The list of prime numbers to check are from the answer choices.
    primes = [17, 383, 1093, 66161, 534851]

    print("We are looking for the smallest prime p in the choices that satisfies the condition: 6^(p-1) ≡ 1 (mod p^2).")
    print("We will test the primes from the list in increasing order.")

    # Iterate through the primes, check the condition, and print the process.
    for p in primes:
        p_squared = p * p
        print(f"\n--- Checking for p = {p} ---")
        
        # We perform the modular exponentiation pow(base, exponent, modulus).
        print(f"Calculating the value of 6^({p - 1}) mod {p}^2, which is 6^({p-1}) mod {p_squared}.")
        
        result = pow(6, p - 1, p_squared)
        
        # Display the full equation with the result.
        print(f"The result is: 6^({p - 1}) mod {p_squared} = {result}")

        if result == 1:
            print(f"Since the result is 1, p = {p} satisfies the condition.")
            print(f"\nThus, the smallest prime in the list for which Z[6^(1/p)] is not the ring of integers is {p}.")
            # Once we find the smallest prime, we stop.
            return
        else:
            print(f"Since the result is {result} (and not 1), p = {p} does not satisfy the condition.")

solve()