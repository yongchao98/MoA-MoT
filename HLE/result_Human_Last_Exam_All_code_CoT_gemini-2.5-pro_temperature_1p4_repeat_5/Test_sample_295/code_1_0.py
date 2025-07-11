def compute_euler_char_mod_k(k):
    """
    Computes the reduced Euler characteristic of the simplicial complex Delta_k, modulo k.

    Args:
        k: A prime number such that k >= 3.
    """
    if not isinstance(k, int) or k < 3:
        print("Error: k must be an integer greater than or equal to 3.")
        return

    # According to the derivation, the reduced Euler characteristic mod k is (k-1)/2.
    # We will show the steps of this simple calculation.
    
    print(f"The given prime is k = {k}.")
    
    # Step 1: Calculate k-1
    k_minus_1 = k - 1
    print(f"The formula requires calculating k-1, which is {k_minus_1}.")

    # Step 2: Divide by 2
    result = k_minus_1 // 2
    print(f"The final result is ({k}-1) / 2 = {k_minus_1} / 2 = {result}.")
    
    # Final answer as per requested format
    print(f"\nFinal Answer: {result}")
    
# You can test this function with a prime k >= 3.
# For example, let's use k = 13.
k_val = 13 
compute_euler_char_mod_k(k_val)
print("\n--- Another example with k=7 ---")
compute_euler_char_mod_k(7)