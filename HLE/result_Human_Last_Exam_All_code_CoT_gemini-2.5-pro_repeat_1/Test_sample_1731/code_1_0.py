import numpy as np

def calculate_and_display_equilibrium_values():
    """
    Derives and displays the equilibrium values for mean energy and entropy
    for the Bose case of light quanta, based on large deviation principles.
    """
    # Define symbolic variables for printing equations
    beta_sym = "β"
    epsilon_sym = "ε"
    k_b_sym = "k"

    # --- Part 1: Symbolic Derivation Display ---
    print("---[ Symbolic Derivation ]---\n")

    # Derivation of Mean Energy
    print("1. Equilibrium Mean Energy (<u>) per Particle:")
    print("---------------------------------------------")
    print(f"The partition function Z for a single mode is the sum of Boltzmann factors over all energy levels ε_j = j*ε (for j=0, 1, 2,...):")
    print(f"Z = Σ_j exp(-{beta_sym}*j*{epsilon_sym})")
    print(f"This is a geometric series which sums to:")
    print(f"Z = 1 / (1 - exp(-{beta_sym}*{epsilon_sym}))\n")

    print(f"The mean energy <u> is derived from the partition function via <u> = -∂(ln(Z)) / ∂{beta_sym}.")
    print(f"This yields the final equation for mean energy (Planck's formula):")
    print(f"<u> = {epsilon_sym} / (exp({beta_sym} * {epsilon_sym}) - 1)\n")

    # Derivation of Entropy
    print("2. Equilibrium Entropy (s) per Particle:")
    print("------------------------------------------")
    print("The entropy per particle 's' is given by the formula s = k * (ln(Z) + β*<u>).")
    print("Substituting the expressions for Z and <u>, we get the final equation:")
    term1_s = f"ln(1 / (1 - exp(-{beta_sym}*{epsilon_sym})))"
    term2_s = f"({beta_sym}*{epsilon_sym}) / (exp({beta_sym}*{epsilon_sym}) - 1)"
    print(f"s = {k_b_sym} * [ {term1_s} + {term2_s} ]")

    print("\n" + "="*50 + "\n")

    # --- Part 2: Numerical Example ---
    print("---[ Numerical Example ]---\n")

    # Use normalized units for demonstration
    beta = 1.0      # Represents inverse temperature, β = 1/(k*T)
    epsilon = 1.0   # Represents quantum of energy, ε = h*ν
    k_b = 1.0       # Boltzmann constant in these units

    print(f"Using example values: {beta_sym} = {beta}, {epsilon_sym} = {epsilon}, {k_b_sym} = {k_b}\n")

    # Calculate and display Mean Energy step-by-step
    print("1. Calculating Mean Energy:")
    mean_energy_numerator = epsilon
    exp_val = np.exp(beta * epsilon)
    mean_energy_denominator = exp_val - 1
    mean_energy_val = mean_energy_numerator / mean_energy_denominator
    
    print(f"<u> = {mean_energy_numerator} / (exp({beta} * {epsilon}) - 1)")
    print(f"<u> = {mean_energy_numerator} / ({exp_val:.4f} - 1)")
    print(f"<u> = {mean_energy_numerator} / {mean_energy_denominator:.4f}")
    print(f"<u> = {mean_energy_val:.4f}\n")

    # Calculate and display Entropy step-by-step
    print("2. Calculating Entropy:")
    Z_val = 1 / (1 - np.exp(-beta * epsilon))
    ln_Z_val = np.log(Z_val)
    beta_u_val = beta * mean_energy_val
    entropy_val = k_b * (ln_Z_val + beta_u_val)
    
    print(f"s = {k_b} * [ ln(1 / (1 - exp(-{beta}*{epsilon}))) + {beta} * {mean_energy_val:.4f} ]")
    print(f"s = {k_b} * [ ln({Z_val:.4f}) + {beta_u_val:.4f} ]")
    print(f"s = {k_b} * [ {ln_Z_val:.4f} + {beta_u_val:.4f} ]")
    print(f"s = {k_b} * [ {ln_Z_val + beta_u_val:.4f} ]")
    print(f"s = {entropy_val:.4f}")

if __name__ == "__main__":
    calculate_and_display_equilibrium_values()