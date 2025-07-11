import numpy as np

def calculate_critical_correlation():
    """
    Calculates the critical correlation between two input populations (v and s)
    required to balance synaptic potentiation and depression in a feedforward network.
    
    The function uses a derived formula based on mean-field approximations.
    """
    
    # --- Parameters of the network and simulation ---
    # N_v: Number of neurons in input layer v
    N_v = 100
    
    # N_s: Number of neurons in input layer s
    N_s = 80
    
    # mu: Average firing rate of input neurons (for both v and s)
    mu = 0.1 # Example value (e.g., in units of spikes per time step)
    
    # theta: Heterosynaptic offset constant from the learning rule
    theta = 0.15 # Example value
    
    print("This script calculates the critical input correlation 'ρ_c' to balance plasticity.")
    print("The formula is derived from the balance of expected weight changes under a mean-field approximation.")
    print("\nGiven Parameters:")
    print(f"  Number of neurons in v (N_v): {N_v}")
    print(f"  Number of neurons in s (N_s): {N_s}")
    print(f"  Average input rate (μ): {mu}")
    print(f"  Heterosynaptic offset (θ): {theta}")
    print("-" * 30)

    # --- Calculation ---
    # The formula for the critical correlation coefficient (rho_c) is:
    # rho_c = ( (N_v + N_s)^2 / (2 * N_v * N_s) ) * (theta - mu)

    # Calculate the numerator and denominator for clarity
    numerator = (N_v + N_s)**2 * (theta - mu)
    denominator = 2 * N_v * N_s

    # Calculate the critical correlation
    rho_c = numerator / denominator

    # --- Output the result ---
    print("\nCalculation Steps:")
    # Using 'f-strings' to embed the numerical values directly into the equation string
    print(f"1. Calculate the critical covariance 'C':")
    # For Poisson inputs, variance = mean = mu. The correlation coeff rho = C / (std(v)*std(s)) = C / mu.
    # So, C = rho * mu.
    # The formula derived is for rho directly, assuming rho = C / mu.
    # The underlying formula for Covariance C is: C = mu * rho_c
    covariance = mu * rho_c
    print(f"   C = (({N_v} + {N_s})^2 / (2 * {N_v} * {N_s})) * {mu} * ({theta} - {mu})")
    print(f"   C = ({ (N_v + N_s)**2 }) / ({2 * N_v * N_s}) * {mu * (theta - mu)}")
    print(f"   C = {covariance:.5f}")


    print("\n2. Calculate the critical correlation coefficient 'ρ_c':")
    print(f"   ρ_c = C / μ")
    print(f"   ρ_c = (({N_v} + {N_s})^2 / (2 * {N_v} * {N_s})) * ({theta} - {mu})")
    print(f"   ρ_c = ({ (N_v + N_s)**2 }) / ({denominator}) * ({theta - mu})")
    print(f"   ρ_c = { (N_v + N_s)**2 / denominator:.4f} * {theta - mu}")
    
    print("\n--- Final Result ---")
    print(f"The critical amount of correlation required is ρ_c = {rho_c:.5f}")


# Execute the function
if __name__ == "__main__":
    calculate_critical_correlation()
