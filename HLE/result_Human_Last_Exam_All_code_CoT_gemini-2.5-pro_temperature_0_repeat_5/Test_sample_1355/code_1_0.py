import numpy as np
from scipy.special import gamma

def calculate_conductance_moment(n):
    """
    Calculates the n-th statistical moment of the dimensionless conductance <g^n>
    for a disordered Majorana wire at the critical point.
    
    The analytical formula is: <g^n> = Gamma(n - 0.5) * Gamma(0.5) / (2 * pi * Gamma(n))
    """
    if n <= 0.5:
        raise ValueError("n must be greater than 0.5 for the Gamma function argument to be positive.")
    
    # Gamma(0.5) is equal to sqrt(pi)
    gamma_half = np.sqrt(np.pi)
    
    # Calculate the numerator and denominator of the formula
    numerator = gamma(n - 0.5) * gamma_half
    denominator = 2 * np.pi * gamma(n)
    
    return numerator / denominator

# --- Main Calculation ---

# 1. Calculate the average conductance <g> (the 1st moment, n=1)
n_avg = 1
average_conductance = calculate_conductance_moment(n_avg)

# 2. Calculate the fourth moment of conductance <g^4> (n=4)
n_fourth = 4
fourth_moment_conductance = calculate_conductance_moment(n_fourth)

# 3. Calculate the final ratio
ratio = fourth_moment_conductance / average_conductance

# --- Output the Results ---
# The problem asks to output each number in the final equation.
# The exact fractional values are <g> = 1/2 and <g^4> = 5/32.
# The ratio is (5/32) / (1/2) = 5/16.

print("This script calculates the ratio between the fourth and first statistical moments of conductance.")
print("-" * 80)
print(f"The average conductance <g> (for n=1) is: {average_conductance:.5f}")
print(f"The fourth moment of conductance <g^4> (for n=4) is: {fourth_moment_conductance:.5f}")
print("-" * 80)
print(f"The ratio <g^4> / <g> is calculated as:")
print(f"{fourth_moment_conductance:.5f} / {average_conductance:.5f} = {ratio:.5f}")
print(f"The exact fractional value of the ratio is 5/16.")

<<<5/16>>>