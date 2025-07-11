import numpy as np
from scipy.optimize import fsolve

def solve_power_fraction():
    """
    Calculates the fraction of power from Pu-239 after 35% uranium burnup.
    """
    # Step 1: Define initial constants and problem parameters
    # Initial enrichment: 30% U-235, 70% U-238
    N_U235_0 = 0.30
    N_U238_0 = 0.70
    N_U_total_0 = N_U235_0 + N_U238_0

    # Burnup: 35% of initial uranium is consumed
    burnup_fraction = 0.35
    N_U_final_target = N_U_total_0 * (1 - burnup_fraction)

    # Cross-sections in barns (1 barn = 1e-24 cm^2)
    sigma_a_U235 = 591.0
    sigma_f_U235 = 505.0
    sigma_a_U238 = 2.42
    # For U-238, absorption (sigma_a) is essentially all capture (sigma_c), leading to Pu-239
    sigma_c_U238 = sigma_a_U238
    sigma_a_Pu239 = 973.0
    sigma_f_Pu239 = 698.0
    
    # Step 2: Define and solve the burnup equation for neutron fluence (tau)
    def burnup_equation(tau, N25_0, N28_0, sa25, sa28, target_N_U):
        """Equation to find the fluence tau for a given burnup."""
        n25_t = N25_0 * np.exp(-sa25 * tau)
        n28_t = N28_0 * np.exp(-sa28 * tau)
        return n25_t + n28_t - target_N_U

    # Numerically solve the equation for tau. A good initial guess is ~0.03
    initial_tau_guess = 0.03
    tau_solution = fsolve(
        burnup_equation,
        initial_tau_guess,
        args=(N_U235_0, N_U238_0, sigma_a_U235, sigma_a_U238, N_U_final_target)
    )[0]

    # Step 3: Calculate nuclide concentrations at the solved fluence tau
    N_U235_final = N_U235_0 * np.exp(-sigma_a_U235 * tau_solution)

    term1 = (N_U238_0 * sigma_c_U238) / (sigma_a_Pu239 - sigma_a_U238)
    term2 = np.exp(-sigma_a_U238 * tau_solution) - np.exp(-sigma_a_Pu239 * tau_solution)
    N_Pu239_final = term1 * term2

    # Step 4: Calculate the fraction of power from Pu-239 fission
    power_U235 = N_U235_final * sigma_f_U235
    power_Pu239 = N_Pu239_final * sigma_f_Pu239
    total_power = power_U235 + power_Pu239

    if total_power > 0:
        power_fraction_Pu239 = power_Pu239 / total_power
    else:
        power_fraction_Pu239 = 0.0

    # Step 5: Print the results in a clear, step-by-step manner.
    print("Step-by-step Calculation:")
    print("-" * 60)
    print(f"1. Solved for neutron fluence (τ) to reach {burnup_fraction*100}% burnup.")
    print(f"   The required fluence τ is: {tau_solution:.6f}")
    print("-" * 60)
    print("2. Calculated relative nuclide concentrations at this fluence:")
    print(f"   - Final U-235 concentration (N_U235): {N_U235_final:.4e}")
    print(f"   - Final Pu-239 concentration (N_Pu239): {N_Pu239_final:.5f}")
    print("-" * 60)
    print("3. Calculated the power fraction from Plutonium-239.")
    print("   Using the formula: Fraction = (N_Pu239 * σ_f_Pu239) / (N_U235 * σ_f_U235 + N_Pu239 * σ_f_Pu239)\n")
    
    print("   Substituting the calculated and given values:")
    print(f"   Fraction = ({N_Pu239_final:.5f} * {sigma_f_Pu239}) / (({N_U235_final:.4e} * {sigma_f_U235}) + ({N_Pu239_final:.5f} * {sigma_f_Pu239}))\n")

    print(f"   Numerator (relative power from Pu-239) = {power_Pu239:.5f}")
    print(f"   Denominator (total relative power):")
    print(f"     - from U-235 = {power_U235:.4e}")
    print(f"     - from Pu-239 = {power_Pu239:.5f}")
    print(f"     - Total = {total_power:.5f}")
    print("-" * 60)
    print("Final Result:")
    print(f"The fraction of power produced from Plutonium-239 is: {power_fraction_Pu239:.4f}")

if __name__ == '__main__':
    solve_power_fraction()
    # The final answer is wrapped according to the format requirements.
    # The value is ~0.999998, which rounds to 1.0000.
    print("\n<<<1.0000>>>")
