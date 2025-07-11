import numpy as np
from scipy.optimize import fsolve

def solve_power_fraction():
    """
    Calculates the fraction of power produced by Pu-239 in a nuclear reactor
    after a specific burnup percentage.
    """
    # --- 1. Define constants and initial conditions ---

    # Initial atom fractions, assuming a normalized total of 1 atom of Uranium.
    N_U235_0 = 0.30
    N_U238_0 = 0.70
    N_total_U_0 = N_U235_0 + N_U238_0

    # Microscopic cross-sections in barns (1 barn = 1e-24 cm^2).
    sigma_a_U235 = 591.0
    sigma_f_U235 = 505.0
    sigma_a_U238 = 2.42  # For U-238, absorption is almost entirely capture.
    sigma_a_Pu239 = 973.0
    sigma_f_Pu239 = 698.0

    # The problem states 35% of the initial uranium has been burned.
    burnup_fraction = 0.35

    # --- 2. Define and solve the burnup equation for neutron fluence (Phi) ---

    def burnup_equation(phi):
        """
        This function represents the burnup condition. We solve for phi where this function is zero.
        Total atoms burned = (Initial atoms - Final atoms) for both U isotopes.
        """
        u235_burned = N_U235_0 * (1 - np.exp(-sigma_a_U235 * phi))
        u238_burned = N_U238_0 * (1 - np.exp(-sigma_a_U238 * phi))
        total_burned = u235_burned + u238_burned
        return total_burned - (burnup_fraction * N_total_U_0)

    # We provide an initial guess for the solver. Based on a quick check,
    # the fluence must be high enough to burn a significant amount of U-238,
    # so we guess a value on the order of 10^-2.
    initial_guess_phi = 0.03
    # Use scipy's fsolve to find the root, which is the fluence Phi.
    fluence_phi = fsolve(burnup_equation, initial_guess_phi)[0]

    # --- 3. Calculate the final atom concentrations at the calculated fluence ---

    # Final concentration of U-235
    N_U235_final = N_U235_0 * np.exp(-sigma_a_U235 * fluence_phi)

    # Final concentration of Pu-239 using the Bateman equation solution
    # for the U-238 -> Pu-239 chain (assuming U-239 and Np-239 decay is instant).
    pre_factor = (sigma_a_U238 * N_U238_0) / (sigma_a_Pu239 - sigma_a_U238)
    exp_term_U238 = np.exp(-sigma_a_U238 * fluence_phi)
    exp_term_Pu239 = np.exp(-sigma_a_Pu239 * fluence_phi)
    N_Pu239_final = pre_factor * (exp_term_U238 - exp_term_Pu239)

    # --- 4. Calculate the power fraction from Pu-239 ---

    # Power is proportional to the fission rate (N * sigma_f).
    power_U235 = N_U235_final * sigma_f_U235
    power_Pu239 = N_Pu239_final * sigma_f_Pu239
    total_power = power_U235 + power_Pu239

    # Avoid division by zero if total_power is somehow zero.
    if total_power == 0:
        power_fraction_pu = 0
    else:
        power_fraction_pu = power_Pu239 / total_power

    # --- 5. Print the final equation with all the numbers ---

    print("The fraction of power from Pu-239 is calculated by dividing the power from Pu-239 by the total power.")
    print("Power is proportional to (Number of Atoms) * (Fission Cross-Section).\n")
    print("Calculated final relative number of atoms:")
    print(f"  N_U235  = {N_U235_final:.5e}")
    print(f"  N_Pu239 = {N_Pu239_final:.5e}\n")

    print("The final equation is:")
    print("Fraction = (N_Pu239 * σ_f_Pu239) / (N_U235 * σ_f_U235 + N_Pu239 * σ_f_Pu239)\n")
    
    print(f"Fraction = ({N_Pu239_final:.5e} * {sigma_f_Pu239}) / (({N_U235_final:.5e} * {sigma_f_U235}) + ({N_Pu239_final:.5e} * {sigma_f_Pu239}))")
    
    print(f"Fraction = {power_Pu239:.5f} / ({power_U235:.5f} + {power_Pu239:.5f})")
    
    print(f"Fraction = {power_Pu239:.5f} / {total_power:.5f}")
    
    print(f"Fraction = {power_fraction_pu:.5f}")

    print(f"\nAfter 35% of the uranium is burned, {power_fraction_pu:.1%} of the power is produced from plutonium-239.")
    
    # Final answer in the specified format
    print(f"<<<{power_fraction_pu:.3f}>>>")

if __name__ == "__main__":
    solve_power_fraction()