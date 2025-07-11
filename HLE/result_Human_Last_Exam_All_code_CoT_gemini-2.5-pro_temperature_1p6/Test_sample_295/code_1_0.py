def compute_euler_characteristic_mod_k(k):
    """
    Computes the reduced Euler characteristic of Delta_k modulo k.

    Args:
      k: A prime number k >= 3.
    """
    if not isinstance(k, int) or k < 3:
        print("Error: k must be an integer prime, k >= 3.")
        return

    # According to the derivation, the result is (k-1)/2.
    # The problem asks to output each number in the final equation.
    # The final equation is: result = (k - 1) / 2
    
    print(f"The problem is to compute chi_hat(Delta_k) mod k for k = {k}.")
    print("The derived formula for the result is (k - 1) / 2.")
    print("\nLet's calculate the numbers in this equation:")
    
    k_val = k
    one = 1
    two = 2
    
    numerator = k_val - one
    result = numerator // two
    
    print(f"k = {k_val}")
    print(f"The numerator is k - {one} = {numerator}")
    print(f"The denominator is {two}")
    print(f"The result is {numerator} / {two} = {result}")
    
    print(f"\nThe value of hat_chi(Delta_{k}) mod {k} is {result}.")

# Example computation for k = 13, which is a prime >= 3.
compute_euler_characteristic_mod_k(13)
