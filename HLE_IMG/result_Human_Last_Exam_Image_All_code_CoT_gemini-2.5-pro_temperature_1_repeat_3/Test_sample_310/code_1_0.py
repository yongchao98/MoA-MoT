import math

def calculate_nn(N, m):
    """
    Calculates the number-average degree of polymerization (Nn) for a linear polymer.

    Args:
        N (int): Initial degree of polymerization.
        m (int): An integer from 1 to 5 defining the degree of destruction.

    Returns:
        float: The calculated Nn value.
    """
    # Degree of destruction 'a'
    a = m / 25.0
    
    # Formula for Nn for linear polymers
    # Nn = N / (1 + a * (N - 1))
    numerator = N
    denominator = 1 + a * (N - 1)
    nn = numerator / denominator
    
    print(f"For the missing simulation (linear polymer, m={m}):")
    print(f"Initial degree of polymerization N = {N}")
    print(f"Degree of destruction a = {m}/25 = {a}")
    print(f"The number-average degree of polymerization is Nn = {N} / (1 + {a} * ({N}-1))")
    print(f"Nn = {numerator} / {denominator}")
    print(f"Nn = {nn:.2f}")
    return nn

# --- Main execution ---
# From the plots, the initial degree of polymerization N is 20.
N_initial = 20

# Our analysis indicates that the plots for linear polymers with m=1,2,3,4 are present,
# while the plot for m=5 is missing.
m_missing = 5

# Calculate Nn for the missing simulation
final_answer = calculate_nn(N_initial, m_missing)

# The final answer in the required format
# print(f"<<<{final_answer:.2f}>>>")