import numpy as np

def calculate_ground_state_occupation(mu, epsilon_0, k_B_T):
    """
    Calculates the occupation number of the ground state based on the
    Bose-Einstein distribution.
    
    The equation is: n_0 = 1 / (exp((epsilon_0 - mu) / k_B_T) - 1)
    """
    
    # The chemical potential mu must be less than the ground state energy epsilon_0
    if mu >= epsilon_0:
        return f"Unphysical: mu ({mu}) must be less than epsilon_0 ({epsilon_0})."

    # Calculate the argument of the exponential
    exponent_arg = (epsilon_0 - mu) / k_B_T
    
    # Calculate the denominator
    denominator = np.exp(exponent_arg) - 1
    
    # To avoid division by zero if the denominator is extremely small
    if denominator < 1e-15:
        occupation_number = float('inf')
    else:
        occupation_number = 1 / denominator

    # Print the numbers used in the final equation
    print(f"For mu = {mu:.4f}:")
    print(f"  Equation: n_0 = 1 / (exp(({epsilon_0} - {mu:.4f}) / {k_B_T}) - 1)")
    print(f"  Result: n_0 = {occupation_number:,.2f} particles")
    print("-" * 40)


if __name__ == "__main__":
    # Define physical constants (in arbitrary units for demonstration)
    epsilon_0 = 0.0  # Ground state energy, set to 0 for convenience
    k_B_T = 1.0      # Thermal energy (k_B * T)

    print(f"Demonstrating Bose-Einstein statistics near condensation.")
    print(f"Ground state energy (epsilon_0) = {epsilon_0}")
    print(f"Thermal energy (k_B*T) = {k_B_T}\n")
    
    # Values of chemical potential (mu) approaching epsilon_0 from below
    mu_values = [-1.0, -0.5, -0.1, -0.01, -0.001, -0.0001, 0.0]

    for mu in mu_values:
        calculate_ground_state_occupation(mu, epsilon_0, k_B_T)
