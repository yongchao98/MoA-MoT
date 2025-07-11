def compute_reduced_euler_char_mod_k(k: int):
    """
    Computes the reduced Euler characteristic of the simplicial complex Delta_k, modulo k.

    Args:
        k: A prime number such that k >= 3.
    """
    # First, we check if k is a prime number >= 3 as per the problem statement.
    if not isinstance(k, int) or k < 3:
        print(f"Error: k must be an integer greater than or equal to 3. Received k={k}.")
        return

    # A simple primality test for demonstration.
    is_prime = True
    if k > 3:
        if k % 2 == 0:
            is_prime = False
        else:
            for i in range(3, int(k**0.5) + 1, 2):
                if k % i == 0:
                    is_prime = False
                    break
    
    if not is_prime:
        print(f"Error: k must be a prime number. {k} is not prime.")
        return

    # The problem asks for the reduced Euler characteristic of Delta_k, modulo k.
    # Our derivation shows that chi_hat(Delta_k) = (k-3)/2 (mod k).
    
    # Since k is an odd prime, k-3 is an even integer.
    # The division (k-3)//2 will result in an integer.
    numerator = k - 3
    denominator = 2
    
    result = numerator // denominator

    # The value (k-3)/2 for k>=3 is always non-negative and less than k,
    # so (k-3)/2 mod k is just (k-3)/2.
    
    print(f"For the prime k = {k}:")
    print(f"The reduced Euler characteristic is computed as (k - 3) / 2.")
    print(f"The final equation is: ( {k} - 3 ) / {denominator} = {result}")
    print(f"So, chi_hat(Delta_{k}) mod {k} is {result}.")


# Example usage for a user-provided k.
# You can change the value of k to any prime number >= 3.
# For example, let's use k = 11.
k_value = 11
compute_reduced_euler_char_mod_k(k_value)