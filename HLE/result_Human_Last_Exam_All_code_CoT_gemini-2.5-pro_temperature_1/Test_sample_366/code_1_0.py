import numpy as np
from scipy.optimize import root_scalar

def solve_power_fraction():
    """
    Calculates the fraction of power produced by Pu-239 after 35% uranium burnup.
    """
    # 1. Define constants based on the problem statement
    
    # Initial fractions and relative numbers (basis of 100 atoms)
    f_U235_0 = 0.30
    f_U238_0 = 0.70
    N_U235_0 = 30.0
    N_U238_0 = 70.0
    N_U_total_0 = N_U235_0 + N_U238_0

    # Burnup fraction
    burnup_fraction = 0.35
    uranium_remaining_fraction = 1.0 - burnup_fraction

    # Cross-sections in barns
    sigma_a_U235 = 591.0   # Absorption cross-section for U-235
    sigma_a_U238 = 2.42    # Absorption cross-section for U-238 (assumed to be capture)
    sigma_a_Pu239 = 973.0  # Absorption cross-section for Pu-239
    sigma_f_U235 = 505.0   # Fission cross-section for U-235
    sigma_f_Pu239 = 698.0  # Fission cross-section for Pu-239

    # 2. Define the function to find the required neutron fluence (phi)
    def find_fluence_equation(fluence):
        """
        This function represents the burnup condition. It should equal zero at the correct fluence.
        The equation is: (N_U235(phi) + N_U238(phi)) / N_U_total_0 - (1 - burnup) = 0
        """
        # Calculate remaining U-235 and U-238 fractions
        rem_U235_frac = f_U235_0 * np.exp(-sigma_a_U235 * fluence)
        rem_U238_frac = f_U238_0 * np.exp(-sigma_a_U238 * fluence)
        
        return rem_U235_frac + rem_U238_frac - uranium_remaining_fraction

    # 3. Solve for the neutron fluence
    # We need a bracket [a, b] where the function value changes sign.
    # f(0) = 0.3 + 0.7 - 0.65 = 0.35 > 0
    # f(0.1) = 0.3*e^(-59.1) + 0.7*e^(-0.242) - 0.65 ~ -0.1 < 0
    # The solution lies between 0 and 0.1.
    try:
        sol = root_scalar(find_fluence_equation, bracket=[0, 0.1], method='brentq')
        fluence = sol.root
    except (ValueError, RuntimeError) as e:
        print(f"Error: Could not solve for fluence. {e}")
        return

    # 4. Calculate final atom numbers using the determined fluence
    N_U235_final = N_U235_0 * np.exp(-sigma_a_U235 * fluence)
    
    # The amount of Pu-239 is calculated from the transmutation of U-238
    # N_Pu239(t) = [N_U238_0 * lambda_U238 / (lambda_Pu239 - lambda_U238)] * (e^(-lambda_U238*t) - e^(-lambda_Pu239*t))
    # where lambda = sigma_a * phi. So we use lambda/phi = sigma_a
    prefactor_pu = (N_U238_0 * sigma_a_U238) / (sigma_a_Pu239 - sigma_a_U238)
    exp_term_pu = np.exp(-sigma_a_U238 * fluence) - np.exp(-sigma_a_Pu239 * fluence)
    N_Pu239_final = prefactor_pu * exp_term_pu

    # 5. Calculate the relative power components
    # Power is proportional to N * sigma_f
    power_U235 = N_U235_final * sigma_f_U235
    power_Pu239 = N_Pu239_final * sigma_f_Pu239

    # 6. Calculate the final fraction of power from Plutonium
    total_power = power_U235 + power_Pu239
    if total_power == 0:
        fraction_Pu_power = 0
    else:
        fraction_Pu_power = power_Pu239 / total_power

    # 7. Print the results clearly
    print("The calculation for the fraction of power from Plutonium-239 is based on the formula:")
    print("Fraction = (N_Pu239 * sigma_f_Pu239) / (N_U235 * sigma_f_U235 + N_Pu239 * sigma_f_Pu239)\n")
    print(f"After 35% uranium burnup, the required neutron fluence is found to be ~{fluence:.5f} barns^-1.\n")
    print("At this point, the relative amounts of key isotopes and their power contributions are:")
    print(f"Final relative number of U-235 atoms: {N_U235_final:.4g}")
    print(f"Final relative number of Pu-239 atoms: {N_Pu239_final:.4g}\n")

    print(f"The relative power component from U-235 (N_U235 * sigma_f_U235) is: {power_U235:.4g}")
    print(f"The relative power component from Pu-239 (N_Pu239 * sigma_f_Pu239) is: {power_Pu239:.4g}\n")
    
    print("The final equation with these values is:")
    print(f"Fraction = {power_Pu239:.4g} / ({power_U235:.4g} + {power_Pu239:.4g})\n")
    
    print(f"The fraction of power produced from Plutonium-239 is: {fraction_Pu_power:.6f}")
    
    # Return the final answer in the required format
    return fraction_Pu_power

if __name__ == '__main__':
    final_answer = solve_power_fraction()
    # The final answer is wrapped in <<<>>> as requested.
    print(f"\n<<<{final_answer}>>>")
