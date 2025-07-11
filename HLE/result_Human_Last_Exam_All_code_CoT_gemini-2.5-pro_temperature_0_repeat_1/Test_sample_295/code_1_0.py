def compute_euler_char_mod_k(k: int):
    """
    Computes the reduced Euler characteristic of the complex Delta_k modulo k.

    Args:
        k: A prime number k >= 3.
    """
    if k < 3 or k % 2 == 0 and k != 2:
        # A simple check, not a full primality test.
        print(f"Warning: k={k} is expected to be a prime >= 3.")

    # The reduced Euler characteristic of Delta_k modulo k is (k-1)/2.
    
    # Step 1: Calculate the numerator of the expression.
    numerator = k - 1

    # Step 2: The denominator is 2.
    denominator = 2

    # Step 3: Compute the result.
    # Since k is an odd prime, k-1 is even, so the division results in an integer.
    result = numerator // denominator

    # Print the final equation and the result.
    print(f"For the prime k = {k}:")
    print(f"The value of the reduced Euler characteristic of Delta_k modulo k is given by the formula (k - 1) / 2.")
    print(f"The calculation is:")
    print(f"({k} - 1) / {denominator} = {numerator} / {denominator} = {result}")
    print(f"Final Answer: The value of chi_hat(Delta_{k}) mod {k} is {result}.")

# Example usage with a prime k >= 3.
# Let's use k = 17 as an example.
compute_euler_char_mod_k(17)