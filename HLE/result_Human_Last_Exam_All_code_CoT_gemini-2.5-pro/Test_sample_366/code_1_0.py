import numpy as np
from scipy.optimize import fsolve

def solve_power_fraction():
    """
    Calculates the fraction of power from Pu-239 in a nuclear reactor
    after a specified level of uranium burnup.
    """
    # --- Step 1: Define Constants from the problem statement ---

    # Initial fuel composition (relative fractions)
    f_U235_0 = 0.30
    f_U238_0 = 0.70

    # Total uranium burnup fraction specified
    burnup_fraction = 0.35

    # Nuclear cross-sections in units of barns
    sigma_a_U235 = 591.0   # Absorption cross-section for U-235
    sigma_f_U235 = 505.0   # Fission cross-section for U-235
    sigma_a_U238 = 2.42    # Absorption cross-section for U-238
    sigma_a_Pu239 = 973.0  # Absorption cross-section for Pu-239
    sigma_f_Pu239 = 698.0  # Fission cross-section for Pu-239

    # --- Step 2: Define and solve the burnup equation for neutron fluence (Theta) ---

    def burnup_equation(Theta):
        """
        This function represents the total uranium burnup. We need to find the
        fluence Theta where this function equals the target burnup_fraction.
        The equation is: f_U235_0*(1-exp(-sa*T)) + f_U238_0*(1-exp(-sa*T)) - burnup = 0
        """
        # fsolve passes an array, so we extract the first element
        Theta = Theta[0]
        # Consumed U-235
        consumed_U235 = f_U235_0 * (1 - np.exp(-sigma_a_U235 * Theta))
        # Consumed U-238
        consumed_U238 = f_U238_0 * (1 - np.exp(-sigma_a_U238 * Theta))
        
        return consumed_U235 + consumed_U238 - burnup_fraction

    # Provide an initial guess for the solver. The high burnup suggests
    # a significant fluence is needed.
    initial_guess_Theta = [0.03]
    # Use fsolve to find the root, which is the required fluence Theta
    fluence_Theta = fsolve(burnup_equation, initial_guess_Theta)[0]

    # --- Step 3: Calculate the final nuclide concentrations at this fluence ---

    # Concentration of U-235 remaining
    N_U235_final = f_U235_0 * np.exp(-sigma_a_U235 * fluence_Theta)

    # Concentration of Pu-239 present. This accounts for both its production
    # from U-238 and its own consumption (burnup).
    prefactor = (f_U238_0 * sigma_a_U238) / (sigma_a_Pu239 - sigma_a_U238)
    exp_term_diff = np.exp(-sigma_a_U238 * fluence_Theta) - np.exp(-sigma_a_Pu239 * fluence_Theta)
    N_Pu239_final = prefactor * exp_term_diff

    # --- Step 4: Calculate the power fraction ---

    # Relative power is proportional to N * sigma_f
    power_from_U235 = N_U235_final * sigma_f_U235
    power_from_Pu239 = N_Pu239_final * sigma_f_Pu239
    total_power = power_from_U235 + power_from_Pu239
    
    if total_power == 0:
        fraction_from_Pu239 = 0
    else:
        fraction_from_Pu239 = power_from_Pu239 / total_power

    # --- Step 5: Print the results, showing the final equation ---
    
    print("This script calculates the power fraction from Plutonium-239.")
    print("-" * 60)
    print(f"To achieve a {burnup_fraction:.0%} burnup of total uranium, the required")
    print(f"neutron fluence (Theta) is calculated to be: {fluence_Theta:.6f} [barns^-1]")
    print("-" * 60)

    print("The final equation for the power fraction from Plutonium-239 is:")
    print("Fraction = Power_Pu239 / (Power_U235 + Power_Pu239)")
    print("\nWhere Power is proportional to (Concentration * Fission Cross Section):\n")
    
    print(f"1. Calculating Relative Power from U-235:")
    print(f"   Power_U235 = N_U235 * sigma_f_U235")
    print(f"              = {N_U235_final:.4e} * {sigma_f_U235} = {power_from_U235:.4e}")
    
    print(f"\n2. Calculating Relative Power from Pu-239:")
    print(f"   Power_Pu239 = N_Pu239 * sigma_f_Pu239")
    print(f"              = {N_Pu239_final:.6f} * {sigma_f_Pu239} = {power_from_Pu239:.4f}")

    print("\n3. Substituting these values into the fraction equation:")
    print(f"   Fraction = {power_from_Pu239:.4f} / ({power_from_U235:.4e} + {power_from_Pu239:.4f})")
    print(f"   Fraction = {power_from_Pu239:.4f} / {total_power:.4f}")
    
    print("-" * 60)
    print(f"Final Answer: The fraction of power produced from Plutonium-239 is {fraction_from_Pu239:.6f}")
    print("-" * 60)

if __name__ == "__main__":
    solve_power_fraction()