import numpy as np
from scipy.optimize import fsolve

def solve_power_fraction():
    """
    Calculates the fraction of power from Pu-239 after 35% uranium burnup.
    """
    # Step 1: Define initial conditions and constants from the problem description
    # Initial normalized number densities (enrichment is 30%)
    N_U235_initial = 0.30
    N_U238_initial = 0.70
    N_U_initial_total = N_U235_initial + N_U238_initial

    # Microscopic cross-sections in barns
    sa_U235 = 591.0   # absorption for U-235
    sa_U238 = 2.42    # absorption for U-238
    sf_U235 = 505.0   # fission for U-235
    sa_Pu239 = 973.0  # absorption for Pu-239
    sf_Pu239 = 698.0  # fission for Pu-239

    # Burnup condition: 35% of uranium is burned, so 65% remains
    final_U_fraction = 1.0 - 0.35

    # Step 2: Determine the required neutron fluence (tau)
    # We need to solve the equation: N_U235(tau) + N_U238(tau) = 0.65 * N_U_initial_total
    # which is N_U235_initial * exp(-sa_U235 * tau) + N_U238_initial * exp(-sa_U238 * tau) = final_U_fraction
    def find_fluence_equation(tau):
        return N_U235_initial * np.exp(-sa_U235 * tau) + N_U238_initial * np.exp(-sa_U238 * tau) - final_U_fraction

    # Solve the equation numerically for tau
    # An initial guess can be made by noting sa_U235 is large, so the first term goes to zero quickly
    initial_guess_tau = -np.log(final_U_fraction / N_U238_initial) / sa_U238
    fluence_tau = fsolve(find_fluence_equation, initial_guess_tau)[0]

    # Step 3: Calculate final isotope concentrations using the determined fluence
    # Final concentration of U-235
    N_U235_final = N_U235_initial * np.exp(-sa_U235 * fluence_tau)

    # Final concentration of Pu-239. It's produced from U-238 and consumed.
    # The solution to the differential equation dN_Pu/dt = R_prod - R_abs is:
    # N_Pu239(tau) = (sa_U238 * N_U238_initial / (sa_Pu239 - sa_U238)) * (exp(-sa_U238 * tau) - exp(-sa_Pu239 * tau))
    term1 = (sa_U238 * N_U238_initial) / (sa_Pu239 - sa_U238)
    term2 = np.exp(-sa_U238 * fluence_tau) - np.exp(-sa_Pu239 * fluence_tau)
    N_Pu239_final = term1 * term2

    # Step 4: Calculate the power contributions and the final fraction
    # Power is proportional to N * sigma_fission
    power_U235 = N_U235_final * sf_U235
    power_Pu239 = N_Pu239_final * sf_Pu239
    
    total_power = power_U235 + power_Pu239
    
    # Avoid division by zero if total power is somehow zero
    if total_power == 0:
        fraction_Pu = 0
    else:
        fraction_Pu = power_Pu239 / total_power
        
    print(f"To achieve 35% uranium burnup, a neutron fluence (τ) of {fluence_tau:.5f} barns⁻¹ is required.\n")
    print("At this fluence, the final relative concentrations are:")
    print(f"  N_U235  = {N_U235_final:.4g}")
    print(f"  N_Pu239 = {N_Pu239_final:.4g}\n")
    
    print("Power produced is proportional to concentration × fission cross-section (N × σ_f).")
    print(f"  Power from U-235  ∝ {N_U235_final:.4g} × {sf_U235} = {power_U235:.4g}")
    print(f"  Power from Pu-239 ∝ {N_Pu239_final:.4g} × {sf_Pu239} = {power_Pu239:.4g}\n")

    print("The fraction of power from Pu-239 is calculated as:")
    print(f"Fraction = Power_Pu239 / (Power_U235 + Power_Pu239)")
    print(f"Fraction = {power_Pu239:.4g} / ({power_U235:.4g} + {power_Pu239:.4g})")
    print(f"Fraction = {fraction_Pu:.5f}")

solve_power_fraction()
<<<0.99998>>>