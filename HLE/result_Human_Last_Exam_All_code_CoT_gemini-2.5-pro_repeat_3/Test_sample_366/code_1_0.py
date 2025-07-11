import numpy as np
from scipy.optimize import fsolve

def solve_power_fraction():
    """
    Calculates the fraction of power produced by Pu-239 in a nuclear reactor
    after a specified level of uranium burnup.
    """
    # Step 1: Define constants and initial conditions
    # Cross-sections in barns (1 barn = 1e-24 cm^2)
    sigma_a_U235 = 591.0
    sigma_f_U235 = 505.0
    sigma_a_U238 = 2.42
    sigma_a_Pu239 = 973.0
    sigma_f_Pu239 = 698.0

    # Initial fuel composition (normalized to 1 total U atom)
    initial_enrichment = 0.30
    N_U235_0 = initial_enrichment
    N_U238_0 = 1.0 - initial_enrichment

    # Target burnup
    burnup_fraction = 0.35

    # Step 2: Define the function to find the required neutron fluence (tau)
    # The function is zero when the remaining uranium equals the post-burnup amount.
    def burnup_equation(tau):
        remaining_U = N_U235_0 * np.exp(-sigma_a_U235 * tau) + N_U238_0 * np.exp(-sigma_a_U238 * tau)
        target_remaining_U = 1.0 - burnup_fraction
        return remaining_U - target_remaining_U

    # Step 3: Solve for the fluence tau_f
    # An initial guess is needed for the numerical solver.
    initial_guess_tau = 0.03
    tau_f = fsolve(burnup_equation, initial_guess_tau)[0]
    
    print(f"Calculated fluence (τ) required for 35% burnup: {tau_f:.5f} [barns^-1]\n")

    # Step 4: Calculate final nuclide concentrations using the solved fluence
    # Final concentration of U-235
    N_U235_f = N_U235_0 * np.exp(-sigma_a_U235 * tau_f)
    
    # Final concentration of Pu-239 using the Bateman equation
    # This accounts for its production from U-238 and its own consumption.
    term1 = (N_U238_0 * sigma_a_U238) / (sigma_a_Pu239 - sigma_a_U238)
    term2 = np.exp(-sigma_a_U238 * tau_f) - np.exp(-sigma_a_Pu239 * tau_f)
    N_Pu239_f = term1 * term2

    # Step 5: Calculate the fission rates (proportional to power)
    power_rate_U235 = N_U235_f * sigma_f_U235
    power_rate_Pu239 = N_Pu239_f * sigma_f_Pu239
    total_power_rate = power_rate_U235 + power_rate_Pu239

    # Calculate the final fraction of power from Pu-239
    fraction_from_Pu239 = power_rate_Pu239 / total_power_rate

    # Step 6: Print the components of the final calculation as requested
    print("--- Final State Concentrations (relative to initial total Uranium) ---")
    print(f"Final concentration of U-235 (N_U235_f): {N_U235_f:.4e}")
    print(f"Final concentration of Pu-239 (N_Pu239_f): {N_Pu239_f:.4e}\n")
    
    print("--- Power Fraction Calculation ---")
    print("Fraction = (N_Pu239_f * σ_f^Pu239) / (N_U235_f * σ_f^U235 + N_Pu239_f * σ_f^Pu239)")
    print(f"Fraction = ({N_Pu239_f:.4e} * {sigma_f_Pu239}) / (({N_U235_f:.4e} * {sigma_f_U235}) + ({N_Pu239_f:.4e} * {sigma_f_Pu239}))")
    print(f"Fraction = {power_rate_Pu239:.4f} / ({power_rate_U235:.4e} + {power_rate_Pu239:.4f})")
    print(f"Fraction = {power_rate_Pu239:.4f} / {total_power_rate:.4f}\n")
    print(f"The fraction of power produced by plutonium-239 is: {fraction_from_Pu239:.4f}")
    
    return fraction_from_Pu239

# Execute the function and store the final answer
final_answer = solve_power_fraction()
print(f"\n<<< {final_answer:.4f} >>>")