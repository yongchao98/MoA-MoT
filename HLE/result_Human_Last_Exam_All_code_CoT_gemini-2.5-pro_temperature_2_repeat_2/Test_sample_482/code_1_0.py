import math

def calculate_critical_correlation(N_v, N_s, mu, theta):
    """
    Calculates the critical amount of correlation 'C' required to balance 
    potentiation and depression in the network.

    Args:
        N_v (int): Number of neurons in the input layer v.
        N_s (int): Number of neurons in the input layer s.
        mu (float): The average rate of activation for input neurons.
        theta (float): The heterosynaptic offset constant.

    Returns:
        float: The calculated critical correlation C = <v_i * s_j>.
    """
    
    # Check for valid inputs
    if N_v <= 0 or N_s <= 0:
        print("Error: N_v and N_s must be positive integers.")
        return float('nan')

    # Numerator of the expression for C:
    # (N_v+N_s)^2*mu*theta - (N_v+N_s)*mu - [N_v(N_v-1) + N_s(N_s-1)]*mu^2
    term1 = ((N_v + N_s)**2) * mu * theta
    term2 = (N_v + N_s) * mu
    term3 = (N_v * (N_v - 1) + N_s * (N_s - 1)) * (mu**2)
    numerator = term1 - term2 - term3
    
    # Denominator of the expression for C
    denominator = 2 * N_v * N_s
    
    # Calculate C
    C = numerator / denominator

    # --- Output ---
    print("This script calculates the critical correlation C = <v_i * s_j> required to balance overall potentiation and depression.")
    print("\nGiven the parameters:")
    print(f"  Number of v neurons (N_v): {N_v}")
    print(f"  Number of s neurons (N_s): {N_s}")
    print(f"  Average input rate (mu):   {mu}")
    print(f"  Heterosynaptic offset (theta): {theta}")
    
    print("\nThe critical correlation is found using the formula:")
    print("C = { (N_v+N_s)^2*μ*θ - (N_v+N_s)*μ - [N_v(N_v-1) + N_s(N_s-1)]*μ^2 } / (2*N_v*N_s)")
    
    print("\nSubstituting the given values into the equation:")
    # Printing each term in the final equation as requested
    print(f"C = ( (({N_v}+{N_s})**2)*{mu}*{theta} - ({N_v}+{N_s})*{mu} - ({N_v}*({N_v}-1) + {N_s}*({N_s}-1))*({mu}**2) ) / (2*{N_v}*{N_s})")
    print(f"C = ( {term1} - {term2} - {term3} ) / {denominator}")
    print(f"C = {numerator} / {denominator}")

    print(f"\nFinal calculated value for the critical correlation:")
    print(f"C = {C}")
    
    # For large, symmetric networks (N_v=N_s=N >> 1), C simplifies to approximately 2*μ*θ - μ^2
    if N_v == N_s and N_v > 100:
        C_approx = 2*mu*theta - mu**2
        print(f"\nFor large, symmetric networks, this approximates to C ≈ 2*μ*θ - μ^2 = {C_approx}")

    return C

# --- Example Usage ---
# We can use some placeholder values to demonstrate the function.
# Let's assume large populations and a low firing rate.
# The exact value of mu is not given, but it is expected to be small.
# The threshold theta defines the turning point for plasticity.

N_v_val = 1000
N_s_val = 1000
mu_val = 0.01      # Example: 1% average activity
theta_val = 0.015  # Example: heterosynaptic threshold

# Execute the calculation and print the results
calculate_critical_correlation(N_v_val, N_s_val, mu_val, theta_val)