import math

def calculate_sum_approximation(n):
    """
    Calculates the approximation for the sum based on the derived formula for a given n.
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: n must be a positive integer.")
        return

    # The derived approximation formula is: n^2/2 + 1/120 + 1/(252*n)
    
    # Calculate each term of the formula
    term1 = n**2 / 2.0
    term2 = 1.0 / 120.0
    term3 = 1.0 / (252.0 * n)
    
    result = term1 + term2 + term3
    
    # Output the equation with the numbers substituted
    print(f"The approximation formula for n = {n} is:")
    print(f"S(n) ≈ ({n}^2 / 2) + (1 / 120) + (1 / (252 * {n}))")
    
    # Output the value of each term
    print("\nBreaking down the calculation:")
    print(f"Value of term 1 ({n}^2 / 2): {term1}")
    print(f"Value of term 2 (1 / 120): {term2}")
    print(f"Value of term 3 (1 / (252 * {n})): {term3}")
    
    # Output the final result
    print(f"\nFinal Result: {result}")

# --- User input ---
# You can change the value of n here.
n_value = 100
# --- End of user input ---

calculate_sum_approximation(n_value)
