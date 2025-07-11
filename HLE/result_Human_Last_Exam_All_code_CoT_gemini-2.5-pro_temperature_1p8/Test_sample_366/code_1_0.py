import numpy as np
from scipy.optimize import brentq

def solve_reactor_problem():
    """
    Solves for the fraction of power produced from Pu-239 after 35% uranium burnup.
    """
    # Step 1: Define constants and initial conditions
    # Initial relative abundances based on 30% enrichment
    N_U235_0 = 0.30
    N_U238_0 = 0.70

    # Uranium burnup fraction
    burnup_fraction = 0.35

    # Microscopic cross-sections in barns
    sigma_a_U235 = 591.0   # Absorption for U-235
    sigma_f_U235 = 505.0   # Fission for U-235
    sigma_a_U238 = 2.42    # Absorption for U-238
    sigma_a_Pu239 = 973.0  # Absorption for Pu-239
    sigma_f_Pu239 = 698.0  # Fission for Pu-239

    # Step 2: Set up and solve the equation for neutron fluence (theta)
    # The total amount of uranium remaining is N_U235(t) + N_U238(t).
    # We need to find the fluence 'theta' where the total uranium has decreased by 35%.
    N_U_total_initial = N_U235_0 + N_U238_0
    N_U_total_final = N_U_total_initial * (1.0 - burnup_fraction)

    def burnup_equation(theta):
        """Equation representing the total uranium amount at a given fluence theta."""
        N_U235_t = N_U235_0 * np.exp(-sigma_a_U235 * theta)
        N_U238_t = N_U238_0 * np.exp(-sigma_a_U238 * theta)
        return (N_U235_t + N_U238_t) - N_U_total_final

    # Find the root of the burnup equation to get the required fluence.
    # The solution for theta must be between 0 and a value where burnup is > 35%.
    # A bracket of [0, 0.1] is sufficient.
    fluence_theta = brentq(burnup_equation, 0, 0.1)

    # Step 3: Calculate the final relative amounts of U-235 and Pu-239
    # Amount of U-235 at the calculated fluence
    N_U235_final = N_U235_0 * np.exp(-sigma_a_U235 * fluence_theta)

    # Amount of Pu-239 at the calculated fluence
    # This is the solution to the differential equation for Pu-239 production and destruction.
    term1 = (N_U238_0 * sigma_a_U238) / (sigma_a_Pu239 - sigma_a_U238)
    term2 = np.exp(-sigma_a_U238 * fluence_theta) - np.exp(-sigma_a_Pu239 * fluence_theta)
    N_Pu239_final = term1 * term2

    # Step 4: Calculate the terms for power contribution
    # Power is proportional to N * sigma_f (relative fission rate)
    power_term_U235 = N_U235_final * sigma_f_U235
    power_term_Pu239 = N_Pu239_final * sigma_f_Pu239
    
    # Step 5: Calculate and print the final power fraction
    power_fraction_Pu239 = power_term_Pu239 / (power_term_U235 + power_term_Pu239)

    print("The fraction of power from Pu-239 is calculated using the formula:")
    print("Fraction = (N_Pu * sigma_f_Pu) / (N_U * sigma_f_U + N_Pu * sigma_f_Pu)")
    print("\nWhere:")
    print(f"  N_Pu (relative amount of Pu-239) = {N_Pu239_final:.5f}")
    print(f"  sigma_f_Pu (fission cross-section of Pu-239) = {sigma_f_Pu239:.1f} barns")
    print(f"  N_U (relative amount of U-235) = {N_U235_final:.5f}")
    print(f"  sigma_f_U (fission cross-section of U-235) = {sigma_f_U235:.1f} barns")
    print("\nPlugging in the numbers:")
    print(f"Fraction = ({N_Pu239_final:.5f} * {sigma_f_Pu239:.1f}) / (({N_U235_final:.5f} * {sigma_f_U235:.1f}) + ({N_Pu239_final:.5f} * {sigma_f_Pu239:.1f}))")
    print(f"Fraction = {power_term_Pu239:.5f} / ({power_term_U235:.5f} + {power_term_Pu239:.5f})")
    
    print("\n--- Final Answer ---")
    print(f"The fraction of power produced from plutonium-239 is: {power_fraction_Pu239:.4f}")
    
    # Return the final numerical answer in the specified format
    print(f"\n<<<{power_fraction_Pu239:.4f}>>>")

# Run the calculation
solve_reactor_problem()