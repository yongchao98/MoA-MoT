import numpy as np
from scipy.optimize import fsolve

def solve_adsorption_model():
    """
    Solves the multi-layer adsorption model using mean-field theory
    and prints the calculation for the average number of layers.
    """
    # 1. Define model parameters based on the problem and assumptions
    # Assumption: Maximum number of layers k_max = 2
    k_max = 2

    # Given parameters (using k_B*T as the unit of energy, so beta=1)
    epsilon_1 = 0.1
    mu = 0.15
    z_l = 4

    # Lateral interaction energy, dependent on k_max
    # ε_l = (0.02)^k_max * k_B*T
    epsilon_l = (0.02)**k_max

    # Assumption: Vertical interaction energy for layers m > 1 is equal to the chemical potential
    # This corresponds to equilibrium with a bulk liquid/condensed phase.
    # ε_m = μ for m >= 2
    epsilon_v = mu
    
    # 2. Define the system of self-consistency equations for fsolve
    # The function F(theta) = theta - G(theta), where G calculates new thetas from old ones.
    # fsolve will find theta such that F(theta) = 0.
    def self_consistency_equations(theta):
        """
        System of equations for self-consistent MFA.
        Args:
            theta (np.ndarray): Array of layer coverages [theta_1, theta_2, ..., theta_k_max]
        Returns:
            np.ndarray: The residuals F(theta) for the solver.
        """
        k = len(theta)
        
        # Calculate the energy E_j for a site being in state j (having j layers)
        # These are dimensionless energies (E / k_B*T)
        energies = np.zeros(k + 1)  # E_0, E_1, ..., E_k
        current_total_energy = 0
        for j in range(1, k + 1):
            eps_j = epsilon_1 if j == 1 else epsilon_v
            mfa_term = z_l * epsilon_l * theta[j-1]
            h_j = eps_j + mfa_term - mu
            current_total_energy += h_j
            energies[j] = current_total_energy

        # Calculate the probabilities P_j of being in state j
        boltzmann_factors = np.exp(-energies)
        z_site = np.sum(boltzmann_factors)
        if z_site == 0:
            # Avoid division by zero, though unlikely here
            probabilities = np.zeros(k + 1)
        else:
            probabilities = boltzmann_factors / z_site
        
        # Calculate the new coverages based on these probabilities
        # theta_m is the probability of having m or more layers
        new_theta = np.zeros(k)
        for m in range(1, k + 1):
            new_theta[m-1] = np.sum(probabilities[m:])
            
        return theta - new_theta

    # 3. Solve the system for the self-consistent coverages
    initial_guess = np.full(k_max, 0.5)
    solution_theta = fsolve(self_consistency_equations, initial_guess)

    # 4. Calculate final quantities using the solved coverages
    final_energies = np.zeros(k_max + 1)
    current_total_energy = 0
    for j in range(1, k_max + 1):
        eps_j = epsilon_1 if j == 1 else epsilon_v
        mfa_term = z_l * epsilon_l * solution_theta[j-1]
        h_j = eps_j + mfa_term - mu
        current_total_energy += h_j
        final_energies[j] = current_total_energy
    
    final_boltzmann_factors = np.exp(-final_energies)
    final_z_site = np.sum(final_boltzmann_factors)
    
    # The average number of layers <n> = sum(j * P_j) = sum(j * exp(-E_j)) / sum(exp(-E_j))
    numerator = np.sum(np.arange(k_max + 1) * final_boltzmann_factors)
    denominator = final_z_site
    avg_layers = numerator / denominator

    # 5. Print the final equation with all numbers included
    print("Based on the assumptions (k_max=2, ε_m=μ for m>1), the average number of layers <n> is calculated.")
    print("The energy of a site with j layers is E_j. The probability is proportional to exp(-E_j/k_B T).")
    print(f"The solved self-consistent layer coverages are: θ₁={solution_theta[0]:.4f}, θ₂={solution_theta[1]:.4f}")
    print(f"The corresponding dimensionless site energies are: E₀/k_BT=0, E₁/k_BT={final_energies[1]:.4f}, E₂/k_BT={final_energies[2]:.4f}")
    
    print("\nThe equation for the average number of layers <n> is:")
    print("<n> = (0*exp(-E₀/k_BT) + 1*exp(-E₁/k_BT) + 2*exp(-E₂/k_BT)) / (exp(-E₀/k_BT) + exp(-E₁/k_BT) + exp(-E₂/k_BT))\n")
    
    print("Plugging in the numbers:")
    numerator_str = f"(0 * {final_boltzmann_factors[0]:.4f} + 1 * {final_boltzmann_factors[1]:.4f} + 2 * {final_boltzmann_factors[2]:.4f})"
    denominator_str = f"({final_boltzmann_factors[0]:.4f} + {final_boltzmann_factors[1]:.4f} + {final_boltzmann_factors[2]:.4f})"
    print(f"<n> = {numerator_str} / {denominator_str}")
    
    print(f"<n> = {numerator:.4f} / {denominator:.4f}")
    print(f"<n> = {avg_layers:.5f}")
    
    # Output the final answer in the required format
    print(f"\n<<<{avg_layers:.5f}>>>")

# Run the solver and print the results
solve_adsorption_model()