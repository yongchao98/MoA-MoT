import numpy as np
from scipy.optimize import root_scalar

def solve_power_fraction():
    """
    Calculates the fraction of power from Pu-239 in a nuclear reactor
    after a specified level of uranium burnup.
    """
    # --- Step 1: Define constants and initial conditions ---
    print("--- Initial Parameters ---")
    # Microscopic cross-sections in barns (1 barn = 1e-24 cm^2)
    s_a_U235 = 591.0   # Absorption cross-section for U-235
    s_f_U235 = 505.0   # Fission cross-section for U-235
    s_a_U238 = 2.42    # Absorption cross-section for U-238
    s_a_Pu239 = 973.0  # Absorption cross-section for Pu-239
    s_f_Pu239 = 698.0  # Fission cross-section for Pu-239

    # Initial atomic fractions in the fuel
    N0_U235_frac = 0.30
    N0_U238_frac = 0.70

    # Total fraction of initial uranium atoms that are consumed ("burned")
    burnup_frac = 0.35
    print(f"Initial U-235 enrichment: {N0_U235_frac*100}%")
    print(f"Total uranium burnup: {burnup_frac*100}%\n")


    # --- Step 2: Set up and solve the burnup equation for neutron fluence (Theta) ---
    print("--- Solving for Neutron Fluence (Θ) ---")
    # The total number of consumed uranium atoms (relative to initial total) is:
    # burnup_frac = N0_U235 * (1 - exp(-s_a_U235*Θ)) + N0_U238 * (1 - exp(-s_a_U238*Θ))
    # We need to find the value of Θ that satisfies this equation.
    def burnup_equation(Theta):
        u235_consumed = N0_U235_frac * (1 - np.exp(-s_a_U235 * Theta))
        u238_consumed = N0_U238_frac * (1 - np.exp(-s_a_U238 * Theta))
        return u235_consumed + u238_consumed - burnup_frac

    # Find the root of the equation using a numerical solver.
    # We bracket the solution between 0 and 0.1 barns^-1.
    try:
        solution = root_scalar(burnup_equation, bracket=[0, 0.1], method='brentq')
        Theta = solution.root
        print(f"Solved Neutron Fluence (Θ) that results in {burnup_frac*100}% burnup: {Theta:.5f} [barns^-1]\n")
    except ValueError:
        print("Error: Could not find a solution for Theta in the given bracket.")
        return


    # --- Step 3: Calculate the final number of atoms for fissile isotopes ---
    print("--- Calculating Final Atom Concentrations ---")
    # The number of remaining U-235 atoms is N(Θ) = N(0) * exp(-σ_a * Θ)
    N_final_U235 = N0_U235_frac * np.exp(-s_a_U235 * Theta)
    print(f"Final relative number of U-235 atoms (N_U235): {N_final_U235:.5f}")

    # The number of Pu-239 atoms is found by solving the differential equation for its production and destruction.
    # The solution (Bateman equation) is:
    # N_Pu239(Θ) = (N0_U238 * σ_a_U238 / (σ_a_Pu239 - σ_a_U238)) * (exp(-σ_a_U238*Θ) - exp(-σ_a_Pu239*Θ))
    coeff = (N0_U238_frac * s_a_U238) / (s_a_Pu239 - s_a_U238)
    exp_term = np.exp(-s_a_U238 * Theta) - np.exp(-s_a_Pu239 * Theta)
    N_final_Pu239 = coeff * exp_term
    print(f"Final relative number of Pu-239 atoms (N_Pu239): {N_final_Pu239:.5f}\n")


    # --- Step 4: Calculate the power fraction from Pu-239 ---
    print("--- Calculating Power Fraction from Plutonium-239 ---")
    # Power is proportional to the fission reaction rate (R = N * σ_f * flux)
    # The flux term cancels out when taking the ratio.
    power_U235 = N_final_U235 * s_f_U235
    power_Pu239 = N_final_Pu239 * s_f_Pu239
    total_power = power_U235 + power_Pu239
    power_frac_Pu239 = power_Pu239 / total_power
    
    # Display the final equation with all the numbers
    print("The fraction is calculated as: (Power_Pu239) / (Power_U235 + Power_Pu239)")
    print("Where Power is proportional to (N_atoms * sigma_fission)\n")
    
    print("Final Equation:")
    print(f"Fraction = (N_Pu239 * σ_f_Pu239) / (N_U235 * σ_f_U235 + N_Pu239 * σ_f_Pu239)")
    print(f"Fraction = ({N_final_Pu239:.5f} * {s_f_Pu239}) / (({N_final_U235:.5f} * {s_f_U235}) + ({N_final_Pu239:.5f} * {s_f_Pu239}))")
    print(f"Fraction = {power_Pu239:.3f} / ({power_U235:.3f} + {power_Pu239:.3f})")
    print(f"Fraction = {power_Pu239:.3f} / {total_power:.3f}\n")

    print(f"Final Answer: The fraction of power being produced from plutonium-239 is {power_frac_Pu239:.4f}")
    
    # Final answer in specified format
    print(f"<<<{power_frac_Pu239:.4f}>>>")

if __name__ == '__main__':
    solve_power_fraction()