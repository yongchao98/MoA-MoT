import math

def calculate_T_norm_for_Jn(n):
    """
    Calculates the 1-norm of the correlation matrix T for the state Jn.
    The formula is valid for odd integers n.

    Args:
        n (int): An odd integer.
    """
    if not isinstance(n, int) or n <= 0 or n % 2 == 0:
        print("Error: Please provide a positive odd integer for n.")
        return

    # Based on the derivation, the 1-norm of the correlation matrix T is 2^(n+1).
    # Step 1: Calculate the exponent
    exponent = n + 1
    
    # Step 2: Calculate the result of 2 to the power of the exponent
    # Using math.pow returns a float, so we cast to int.
    # Alternatively, the ** operator can be used for integer exponentiation.
    result = 2**exponent

    # Output the result in a descriptive equation format.
    print(f"For the given odd n = {n}:")
    print(f"The 1-norm of the correlation matrix T is given by the formula 2^(n+1).")
    print(f"Calculation: 2^({n}+1) = 2^{exponent} = {result}")

# Example usage of the function for an arbitrary odd n, e.g., n=5.
# You can change this value to any positive odd integer.
odd_n = 5
calculate_T_norm_for_Jn(odd_n)

# Another example for n=1
calculate_T_norm_for_Jn(1)

# And for n=3
calculate_T_norm_for_Jn(3)