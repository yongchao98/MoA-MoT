import numpy as np

def calculate_partition_function(mu, beta, M, particle_type):
    """
    Calculates the grand canonical partition function Z for a system with H=0.
    This corresponds to the derived formulas for non-interacting particles.

    Args:
        mu (float): Chemical potential. Must be negative for bosons.
        beta (float): Inverse temperature (1/kT). Must be positive.
        M (int): The number of available single-particle states (assumed degenerate).
        particle_type (str): 'boson' or 'fermion'.
    """
    print(f"--- {particle_type.capitalize()} System ---")
    if particle_type.lower() == 'boson':
        if mu >= 0:
            print("Error: For bosons, the chemical potential (mu) must be negative for the partition function to converge.")
            return

        # For M degenerate states, Z = (1 / (1 - exp(beta*mu)))^M
        term = np.exp(beta * mu)
        z_k = 1 / (1 - term)
        # Using np.power for potential large exponents
        Z = np.power(z_k, M)

        print(f"The formula for a single state is: Z_k = 1 / (1 - exp(βμ))")
        print("For M degenerate states, the total partition function is Z = (Z_k)^M")
        print("\nCalculating with the given values:")
        print(f"μ (chemical potential) = {mu}")
        print(f"β (inverse temperature) = {beta}")
        print(f"M (number of states) = {M}")

        # Outputting each component of the final equation
        print("\nEquation with numbers:")
        print(f"Z = (1 / (1 - exp({beta} * {mu})))**{M}")
        print(f"Z = (1 / (1 - {term:.4f}))**{M}")
        print(f"Z = ({z_k:.4f})**{M}")
        print(f"Final result: Z = {Z}")

    elif particle_type.lower() == 'fermion':
        # For M degenerate states, Z = (1 + exp(beta*mu))^M
        term = np.exp(beta * mu)
        z_k = 1 + term
        # Using np.power for potential large exponents
        Z = np.power(z_k, M)
        
        print(f"The formula for a single state is: Z_k = 1 + exp(βμ)")
        print("For M degenerate states, the total partition function is Z = (Z_k)^M")
        print("\nCalculating with the given values:")
        print(f"μ (chemical potential) = {mu}")
        print(f"β (inverse temperature) = {beta}")
        print(f"M (number of states) = {M}")

        # Outputting each component of the final equation
        print("\nEquation with numbers:")
        print(f"Z = (1 + exp({beta} * {mu}))**{M}")
        print(f"Z = (1 + {term:.4f})**{M}")
        print(f"Z = ({z_k:.4f})**{M}")
        print(f"Final result: Z = {Z}")

    else:
        print(f"Error: Unknown particle type '{particle_type}'. Please use 'boson' or 'fermion'.")

# --- Example Usage ---
# You can change these values to explore different scenarios.
# Note: For bosons, mu must be negative.
example_mu = -0.5
example_beta = 1.0
example_M = 10     # Assuming there are 10 degenerate states

# Calculate for a bosonic system
calculate_partition_function(mu=example_mu, beta=example_beta, M=example_M, particle_type='boson')

print("\n" + "="*50 + "\n")

# Calculate for a fermionic system with the same parameters
calculate_partition_function(mu=example_mu, beta=example_beta, M=example_M, particle_type='fermion')
