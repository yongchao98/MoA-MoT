import numpy as np
from scipy.optimize import root_scalar

def solve_power_fraction():
    """
    Calculates the fraction of power produced by Pu-239 after 35% of the initial
    uranium has been consumed.
    """
    # Given cross-sections in barns
    sigma_a_U235 = 591.0
    sigma_a_U238 = 2.42
    sigma_a_Pu239 = 973.0
    sigma_f_Pu239 = 698.0
    sigma_f_U235 = 505.0

    # Initial fuel composition (normalized to 1 total atom)
    N_U235_0 = 0.30
    N_U238_0 = 0.70
    N_total_0 = N_U235_0 + N_U238_0

    # Target burnup
    burnup_fraction = 0.35
    N_U_final_target = N_total_0 * (1 - burnup_fraction)

    # Define the function to find the root for (fluence * sigma)
    # The function is f(Theta) = N_U_final(Theta) - N_U_final_target = 0
    def burnup_equation(theta):
        n_u235_t = N_U235_0 * np.exp(-sigma_a_U235 * theta)
        n_u238_t = N_U238_0 * np.exp(-sigma_a_U238 * theta)
        return n_u235_t + n_u238_t - N_U_final_target

    # Solve for theta using a numerical root finder.
    # A bracket is chosen based on preliminary analysis showing the root is between 0.03 and 0.04.
    solution = root_scalar(burnup_equation, bracket=[0.03, 0.04])
    theta_final = solution.root

    # Calculate final number densities of U-235 and Pu-239
    N_U235_final = N_U235_0 * np.exp(-sigma_a_U235 * theta_final)

    term1 = (N_U238_0 * sigma_a_U238) / (sigma_a_Pu239 - sigma_a_U238)
    term2 = np.exp(-sigma_a_U238 * theta_final) - np.exp(-sigma_a_Pu239 * theta_final)
    N_Pu239_final = term1 * term2

    # Calculate the power contribution from each nuclide (proportional to N * sigma_f)
    power_U235 = N_U235_final * sigma_f_U235
    power_Pu239 = N_Pu239_final * sigma_f_Pu239

    # Calculate the fraction of power from Pu-239
    total_power = power_U235 + power_Pu239
    fraction_Pu_power = power_Pu239 / total_power

    # Print the final equation with all values
    print("The fraction of power produced by Plutonium-239 is calculated as:")
    print("Fraction = (N_Pu239 * σ_f_Pu239) / (N_U235 * σ_f_U235 + N_Pu239 * σ_f_Pu239)\n")
    print(f"Where at 35% burnup:")
    print(f"  N_U235 (relative abundance) = {N_U235_final:.6e}")
    print(f"  N_Pu239 (relative abundance) = {N_Pu239_final:.6f}")
    print(f"  σ_f_U235 = {sigma_f_U235} barns")
    print(f"  σ_f_Pu239 = {sigma_f_Pu239} barns\n")

    print("Substituting the values into the equation:")
    print(f"Fraction = ({N_Pu239_final:.6f} * {sigma_f_Pu239}) / (({N_U235_final:.6e} * {sigma_f_U235}) + ({N_Pu239_final:.6f} * {sigma_f_Pu239}))")
    print(f"Fraction = {power_Pu239:.6f} / ({power_U235:.6e} + {power_Pu239:.6f})")
    print(f"Fraction = {power_Pu239:.6f} / {total_power:.6f}")
    print(f"Final Answer: {fraction_Pu_power:.6f}")
    print(f"\n<<<${fraction_Pu_power:.4f}$>>>")


solve_power_fraction()