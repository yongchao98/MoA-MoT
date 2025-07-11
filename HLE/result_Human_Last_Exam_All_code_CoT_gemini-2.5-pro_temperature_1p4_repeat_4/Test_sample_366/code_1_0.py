import numpy as np
from scipy.optimize import root_scalar

def solve_power_fraction():
    """
    Calculates the fraction of power from Pu-239 after 35% uranium burnup.
    """
    # Step 1: Define physical constants and initial conditions from the problem.
    # We can use relative initial numbers based on the 30% enrichment.
    N0_U235_frac = 0.30  # Initial fraction of U-235
    N0_U238_frac = 0.70  # Initial fraction of U-238

    # Burnup condition
    burnup_frac = 0.35
    remaining_U_frac = 1.0 - burnup_frac

    # Microscopic cross-sections (in barns)
    sigma_a_U235 = 591.0
    sigma_f_U235 = 505.0
    sigma_a_U238 = 2.42
    sigma_a_Pu239 = 973.0
    sigma_f_Pu239 = 698.0

    # Step 2: Set up the equation for fluence Theta based on the burnup condition.
    # The equation is: N0_U235_frac * exp(-sigma_a_U235*Theta) + N0_U238_frac * exp(-sigma_a_U238*Theta) = remaining_U_frac
    def find_fluence_equation(Theta):
        """Equation that equals zero when Theta corresponds to the target burnup."""
        term_U235 = N0_U235_frac * np.exp(-sigma_a_U235 * Theta)
        term_U238 = N0_U238_frac * np.exp(-sigma_a_U238 * Theta)
        return term_U235 + term_U238 - remaining_U_frac

    # Step 3: Solve for Theta numerically using a root-finding algorithm.
    # We search for a root in the interval [0, 1], which is a reasonable range for Theta.
    solution = root_scalar(find_fluence_equation, bracket=[0, 1])
    Theta = solution.root

    # Step 4: Calculate the final relative number of atoms for U-235 and Pu-239.
    # N_U235(Theta) = N_U235_0 * exp(-sigma_a_U235 * Theta)
    N_U235 = N0_U235_frac * np.exp(-sigma_a_U235 * Theta)

    # N_Pu239(Theta) is calculated from U-238 capture, followed by decays.
    factor = (N0_U238_frac * sigma_a_U238) / (sigma_a_Pu239 - sigma_a_U238)
    exp_term_U238 = np.exp(-sigma_a_U238 * Theta)
    exp_term_Pu239 = np.exp(-sigma_a_Pu239 * Theta)
    N_Pu239 = factor * (exp_term_U238 - exp_term_Pu239)

    # Step 5: Calculate the fission reaction rates (proportional to power).
    power_U235 = N_U235 * sigma_f_U235
    power_Pu239 = N_Pu239 * sigma_f_Pu239

    # Step 6: Compute the fraction of total power coming from Pu-239.
    total_power = power_U235 + power_Pu239
    fraction_from_Pu239 = power_Pu239 / total_power

    # Step 7: Print the results, showing the final equation and its components.
    print("--- Calculation Steps and Results ---")
    print(f"Solved for neutron fluence (Theta): {Theta:.6f} [barns^-1]")
    print("\n--- Final State of the Fuel ---")
    print(f"Relative number of U-235 atoms (N_U235): {N_U235:.4e}")
    print(f"Relative number of Pu-239 atoms (N_Pu239): {N_Pu239:.6f}")
    
    print("\n--- Power Fraction Calculation ---")
    print("The final calculation is: Fraction = (N_Pu239 * sigma_f_Pu239) / (N_U235 * sigma_f_U235 + N_Pu239 * sigma_f_Pu239)")
    print("\nPlugging in the numbers into the equation:")
    print(f"Fraction = ({N_Pu239:.6f} * {sigma_f_Pu239}) / (({N_U235:.4e} * {sigma_f_U235}) + ({N_Pu239:.6f} * {sigma_f_Pu239}))")
    print(f"Fraction = {power_Pu239:.6f} / ({power_U235:.4e} + {power_Pu239:.6f})")
    print(f"Fraction = {power_Pu239:.6f} / {total_power:.6f}")

    print("\n--- Final Answer ---")
    print(f"Fraction of power produced from Plutonium-239: {fraction_from_Pu239:.6f}")

    return fraction_from_Pu239

if __name__ == '__main__':
    final_fraction = solve_power_fraction()
    # The required output format for the answer
    print(f"\n<<<{final_fraction}>>>")
