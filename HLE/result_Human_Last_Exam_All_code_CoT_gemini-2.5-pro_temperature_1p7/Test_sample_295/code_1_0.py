def compute_euler_char_mod_k(k):
    """
    Computes the reduced Euler characteristic of Delta_k modulo k.

    Args:
      k: A prime number k >= 3.
    """
    if not isinstance(k, int) or k < 3:
        print("Error: k must be an integer greater than or equal to 3.")
        return

    # According to the derivation, the reduced Euler characteristic modulo k
    # is given by the formula (k-1)/2.
    # Since k is an odd prime, k-1 is even, so the division is exact.
    
    numerator = k - 1
    denominator = 2
    
    # In Python, // performs integer division.
    result = numerator // denominator

    print(f"For k = {k}, the computation of hat_chi(Delta_k) mod k is as follows:")
    # The final equation shows the numbers used in the calculation.
    print(f"({k} - 1) / {denominator} = {numerator} / {denominator} = {result}")
    
    # Since the result is less than k, result mod k is just the result.
    print(f"\nThe value of hat_chi(Delta_{k}) mod {k} is {result}.")

# Example usage of the function with a specific prime k.
# You can change this value to any other prime >= 3.
prime_k = 13
compute_euler_char_mod_k(prime_k)
