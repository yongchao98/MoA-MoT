import math

def calculate_norm_T(n):
    """
    Calculates the 1-norm of the correlation matrix T for the state J_n for odd n.
    The analytical derivation shows that for any odd n, the 1-norm is given by
    the simple formula: 2^(n+1) - 1.

    This script demonstrates the calculation for a specific value of n.
    """
    if n <= 0 or n % 2 == 0:
        print(f"Error: The provided n={n} is not a positive odd integer.")
        print("The derivation is valid for positive odd integers only.")
        return

    # The formula for the 1-norm of the correlation matrix T.
    # ||T||_1 = 2^(n+1) - 1

    print(f"For the quantum state J_n with odd n = {n}, the 1-norm of the correlation matrix T is calculated as follows:")
    print(f"The general formula is ||T||_1 = 2^(n+1) - 1.")

    # Step 1: Calculate the exponent
    exponent = n + 1
    print(f"First, we calculate the exponent: n + 1 = {n} + 1 = {exponent}")

    # Step 2: Calculate the power of 2
    power_of_2 = 2**exponent
    print(f"Next, we compute the power of 2: 2^{exponent} = {power_of_2}")

    # Step 3: Subtract 1 to get the final result
    result = power_of_2 - 1
    print(f"Finally, we subtract 1: {power_of_2} - 1 = {result}")
    
    print("\n------------------------------------")
    print(f"The final result is: {result}")
    print("------------------------------------")


# We choose a specific odd integer n to demonstrate the calculation.
# Let's use n = 5.
if __name__ == "__main__":
    n_example = 5
    calculate_norm_T(n_example)