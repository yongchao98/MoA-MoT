import numpy as np
from scipy.optimize import fsolve

def solve_power_fraction():
    """
    Calculates the fraction of power produced from Plutonium-239
    after 35% of the total uranium has been consumed.
    """
    # Step 1: Define initial conditions and constants from the problem.
    # Initial fractions are based on 30% enrichment, normalizing total uranium atoms to 1.
    initial_enrichment = 0.30
    burnup_fraction = 0.35

    N_U235_initial = initial_enrichment
    N_U238_initial = 1.0 - initial_enrichment
    
    # Microscopic cross-sections in barns (1 barn = 1e-24 cm^2)
    sigma_a_U235 = 591.0   # Absorption cross-section for U-235
    sigma_f_U235 = 505.0   # Fission cross-section for U-235
    sigma_a_U238 = 2.42    # Absorption cross-section for U-238
    sigma_a_Pu239 = 973.0  # Absorption cross-section for Pu-239
    sigma_f_Pu239 = 698.0  # Fission cross-section for Pu-239

    # Step 2: Define and solve the equation for neutron fluence (tau).
    # The total number of uranium atoms remaining is (1 - burnup_fraction) of the initial amount.
    # N_U235_final + N_U238_final = N_initial * (1 - burnup_fraction)
    # N_U235_initial*exp(-sigma_a_U235*tau) + N_U238_initial*exp(-sigma_a_U238*tau) = 1.0 * (1 - burnup_fraction)
    def fluence_equation(tau):
        remaining_U235 = N_U235_initial * np.exp(-sigma_a_U235 * tau)
        remaining_U238 = N_U238_initial * np.exp(-sigma_a_U238 * tau)
        return remaining_U235 + remaining_U238 - (1.0 - burnup_fraction)

    # Numerically solve for tau. An initial guess of 0.03 is reasonable.
    initial_guess_tau = 0.03
    fluence_tau = fsolve(fluence_equation, initial_guess_tau)[0]

    print(f"Calculated neutron fluence (τ): {fluence_tau:.6f} [1/barns]")
    print("-" * 30)

    # Step 3: Calculate the number of U-235 and Pu-239 atoms at this fluence.
    N_U235_final = N_U235_initial * np.exp(-sigma_a_U235 * fluence_tau)
    
    # The concentration of Pu-239 is found by solving the Bateman equation for the U-238 -> Pu-239 chain.
    # d(N_Pu239)/d(tau) = N_U238(tau)*sigma_a_U238 - N_Pu239(tau)*sigma_a_Pu239
    # We assume instantaneous decay of intermediate products (U-239, Np-239).
    exp_U238 = np.exp(-sigma_a_U238 * fluence_tau)
    exp_Pu239 = np.exp(-sigma_a_Pu239 * fluence_tau)
    
    term1 = (N_U238_initial * sigma_a_U238) / (sigma_a_Pu239 - sigma_a_U238)
    N_Pu239_final = term1 * (exp_U238 - exp_Pu239)
    
    print(f"Final relative number of U-235 atoms: {N_U235_final:.4e}")
    print(f"Final relative number of Pu-239 atoms: {N_Pu239_final:.6f}")
    print("-" * 30)
    
    # Step 4: Calculate the fission reaction rate terms for each nuclide.
    # Power is proportional to N * sigma_f.
    power_term_U235 = N_U235_final * sigma_f_U235
    power_term_Pu239 = N_Pu239_final * sigma_f_Pu239
    
    # Step 5: Calculate and print the final power fraction.
    total_power_term = power_term_U235 + power_term_Pu239
    if total_power_term == 0:
        power_fraction_Pu239 = 0
    else:
        power_fraction_Pu239 = power_term_Pu239 / total_power_term

    print("Fraction of Power from Plutonium Fission Calculation:")
    print(f"Fraction = (N_Pu239 * σ_f_Pu239) / (N_U235 * σ_f_U235 + N_Pu239 * σ_f_Pu239)")
    print(f"Fraction = ({N_Pu239_final:.6f} * {sigma_f_Pu239}) / (({N_U235_final:.4e} * {sigma_f_U235}) + ({N_Pu239_final:.6f} * {sigma_f_Pu239}))")
    print(f"Fraction = {power_term_Pu239:.6f} / ({power_term_U235:.6f} + {power_term_Pu239:.6f})")
    print(f"Fraction = {power_term_Pu239:.6f} / {total_power_term:.6f}")
    print(f"Final Power Fraction from Pu-239: {power_fraction_Pu239:.6f}")

    # Returning the final numeric answer for the <<<>>> tag.
    return power_fraction_Pu239

if __name__ == "__main__":
    final_answer = solve_power_fraction()
    # Format the final answer into the required tag format.
    # Due to extreme burnup, almost all U-235 is gone, so the power from Pu-239 is nearly 100%.
    # We present the answer rounded to 4 decimal places.
    print(f"\n<<<After 35% of the uranium has been burned, the fraction of the power being produced from plutonium-239 is {final_answer:.4f}.>>>")
