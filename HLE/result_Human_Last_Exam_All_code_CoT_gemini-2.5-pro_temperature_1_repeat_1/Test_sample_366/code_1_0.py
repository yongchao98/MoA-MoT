import numpy as np
from scipy.optimize import root_scalar

def solve_power_fraction():
    """
    Calculates the fraction of power from Pu-239 after a specific uranium burnup.
    """
    # --- Given Constants ---
    # Absorption cross-sections (in barns)
    sigma_a_U235 = 591.0
    sigma_a_U238 = 2.42
    sigma_a_Pu239 = 973.0

    # Fission cross-sections (in barns)
    sigma_f_U235 = 505.0
    sigma_f_Pu239 = 698.0

    # --- Initial Conditions ---
    # Initial enrichment: 30% U-235, 70% U-238.
    # We normalize the total initial number of uranium atoms to 1.
    N0_U235 = 0.30
    N0_U238 = 0.70
    N0_U_total = N0_U235 + N0_U238 # This is 1.0

    # Target condition: 35% of initial uranium is burned.
    # This means 65% of the initial uranium remains.
    target_N_U_total = N0_U_total * (1.0 - 0.35)

    # --- Step 1: Find the neutron fluence (theta) for the given burnup ---
    print("Step 1: Finding the neutron fluence (theta) for 35% uranium burnup.")

    # We need to solve for theta in the equation:
    # N_U235(theta) + N_U238(theta) = target_N_U_total
    # N0_U235*exp(-sigma_a_U235*theta) + N0_U238*exp(-sigma_a_U238*theta) = target_N_U_total
    def burnup_equation(theta):
        return (N0_U235 * np.exp(-sigma_a_U235 * theta) +
                N0_U238 * np.exp(-sigma_a_U238 * theta) -
                target_N_U_total)

    # Find the root of the equation. The root must be between 0.01 and 0.1 barns^-1.
    solution = root_scalar(burnup_equation, bracket=[0.01, 0.1], method='brentq')
    theta = solution.root
    print(f"Solved fluence, theta = {theta:.6f} barns^-1")

    # --- Step 2: Calculate nuclide concentrations at this fluence ---
    print("\nStep 2: Calculating nuclide concentrations at this fluence.")

    # Concentration of U-235
    N_U235_t = N0_U235 * np.exp(-sigma_a_U235 * theta)
    print(f"Final relative concentration of U-235 (N_U235) = {N_U235_t:.4e}")

    # Concentration of Pu-239 is calculated using the Bateman equation, assuming
    # the U-239 and Np-239 decay is instantaneous on this timescale.
    # d(N_Pu239)/d(theta) = N_U238*sigma_a_U238 - N_Pu239*sigma_a_Pu239
    prefactor = (N0_U238 * sigma_a_U238) / (sigma_a_Pu239 - sigma_a_U238)
    exp_term = np.exp(-sigma_a_U238 * theta) - np.exp(-sigma_a_Pu239 * theta)
    N_Pu239_t = prefactor * exp_term
    print(f"Final relative concentration of Pu-239 (N_Pu239) = {N_Pu239_t:.4e}")

    # --- Step 3: Calculate the power fraction from Pu-239 ---
    print("\nStep 3: Calculate the fraction of power from Plutonium-239.")
    print("The power produced is proportional to N * sigma_f.")

    # Calculate the terms for the power fraction equation
    power_U235_comp = N_U235_t * sigma_f_U235
    power_Pu239_comp = N_Pu239_t * sigma_f_Pu239

    # Calculate the final fraction
    total_power_comp = power_U235_comp + power_Pu239_comp
    fraction_Pu = power_Pu239_comp / total_power_comp

    # Print the final calculation showing all numbers as requested
    print("\nFinal power fraction equation:")
    print(f"Fraction = (N_Pu239 * sigma_f_Pu239) / (N_U235 * sigma_f_U235 + N_Pu239 * sigma_f_Pu239)")
    print(f"Fraction = ({N_Pu239_t:.6e} * {sigma_f_Pu239}) / (({N_U235_t:.6e} * {sigma_f_U235}) + ({N_Pu239_t:.6e} * {sigma_f_Pu239}))")
    print(f"Fraction = ({power_Pu239_comp:.6f}) / ({power_U235_comp:.6e} + {power_Pu239_comp:.6f})")
    print(f"\nResultant Fraction: {fraction_Pu:.6f}")

    return fraction_Pu

if __name__ == '__main__':
    final_fraction = solve_power_fraction()
    # The final answer is the fraction, which is very close to 1.
    print(f"\nAfter 35% of the uranium has been burned, the fraction of power produced from Pu-239 is approximately {final_fraction*100:.4f}%.")
    # Return the final numerical answer in the required format.
    print(f"\n<<<{final_fraction:.6f}>>>")
