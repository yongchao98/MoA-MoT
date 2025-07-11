import numpy as np
from scipy.optimize import fsolve

def solve_power_fraction():
    """
    Calculates the fraction of power produced from Pu-239 after 35% uranium burnup.
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
    # Based on 30% U-235 enrichment. We'll use a basis of 100 initial total uranium atoms.
    N0_U235 = 30.0
    N0_U238 = 70.0
    N0_U_total = N0_U235 + N0_U238

    # --- Burnup Condition ---
    # After 35% of the initial uranium has been consumed.
    burnup_fraction = 0.35
    N_final_U_total = N0_U_total * (1 - burnup_fraction)

    # --- Step 1: Find the required neutron fluence (theta) ---
    # We solve for theta in the equation:
    # N_final_U_total = N0_U235 * exp(-sigma_a_U235 * theta) + N0_U238 * exp(-sigma_a_U238 * theta)
    # The function below represents the equation set to zero.
    def find_fluence_eq(theta):
        # We need to access the variables from the outer scope
        theta = theta[0] # fsolve passes theta as an array
        return N0_U235 * np.exp(-sigma_a_U235 * theta) + N0_U238 * np.exp(-sigma_a_U238 * theta) - N_final_U_total

    # We make an initial guess for theta. Since U-235 burns much faster than U-238,
    # the final number of U atoms will be dominated by the remaining U-238.
    # So, 65 ≈ 70 * exp(-2.42 * theta). This gives theta ≈ 0.03.
    initial_guess = [0.03]
    fluence_theta = fsolve(find_fluence_eq, initial_guess)[0]

    # --- Step 2: Calculate the number of atoms at this fluence ---
    # Number of U-235 atoms remaining
    N_final_U235 = N0_U235 * np.exp(-sigma_a_U235 * fluence_theta)

    # Number of Pu-239 atoms created and remaining.
    # This is the solution to the Bateman equation for Pu-239 production and decay.
    # We assume the capture cross-section for U-238 (sigma_c_U238) is equal to its
    # absorption cross-section, as fission is negligible for U-238 at these energies.
    sigma_c_U238 = sigma_a_U238
    
    # Pre-calculating terms for clarity
    coeff = (sigma_c_U238 * N0_U238) / (sigma_a_Pu239 - sigma_c_U238)
    exp_term_U238 = np.exp(-sigma_c_U238 * fluence_theta)
    exp_term_Pu239 = np.exp(-sigma_a_Pu239 * fluence_theta)
    
    N_final_Pu239 = coeff * (exp_term_U238 - exp_term_Pu239)

    # --- Step 3: Calculate the components of power from each fissile isotope ---
    # Power is proportional to the macroscopic fission rate (N * sigma_f * flux).
    # The flux term is common to both and will cancel out in the fraction.
    power_U235 = N_final_U235 * sigma_f_U235
    power_Pu239 = N_final_Pu239 * sigma_f_Pu239
    
    # --- Step 4: Calculate the fraction of power from Pu-239 ---
    total_fission_power = power_U235 + power_Pu239
    
    if total_fission_power == 0:
        fraction_Pu_power = 0.0
    else:
        fraction_Pu_power = power_Pu239 / total_fission_power

    # --- Step 5: Print the detailed calculation and result ---
    print("Calculation for the fraction of power from Plutonium-239:\n")
    print("Fraction = Power_Pu239 / (Power_U235 + Power_Pu239)")
    print("         = (N_Pu239 * σ_f_Pu239) / ( (N_U235 * σ_f_U235) + (N_Pu239 * σ_f_Pu239) )")
    print("-" * 70)
    print(f"Calculated number of atoms at {burnup_fraction*100}% burnup:")
    print(f"  N_U235  = {N_final_U235:.4e} atoms")
    print(f"  N_Pu239 = {N_final_Pu239:.4f} atoms\n")
    print("Substituting the values into the equation:")
    print(f"Fraction = ({N_final_Pu239:.4f} * {sigma_f_Pu239}) / ( ({N_final_U235:.4e} * {sigma_f_U235}) + ({N_final_Pu239:.4f} * {sigma_f_Pu239}) )")
    print(f"         = {power_Pu239:.4f} / ( {power_U235:.4e} + {power_Pu239:.4f} )")
    print(f"         = {power_Pu239:.4f} / {total_fission_power:.4f}")
    print(f"\nFinal Answer = {fraction_Pu_power:.6f}")
    
    # Returning final answer for the wrapper
    return fraction_Pu_power

# Execute the function and print the final answer in the required format
final_fraction = solve_power_fraction()
print(f"\n<<<The fraction of power from plutonium-239 is {final_fraction:.6f}>>>")
