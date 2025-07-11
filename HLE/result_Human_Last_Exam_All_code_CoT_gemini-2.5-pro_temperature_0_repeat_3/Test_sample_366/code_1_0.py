import math

def solve_power_fraction():
    """
    Calculates the fraction of power produced from plutonium-239
    in a nuclear reactor after 35% of the initial uranium has been burned.
    """
    # --- Step 0: Define constants and initial conditions ---
    # Initial number fractions (normalized to 1 total U atom)
    N_U235_0 = 0.30
    N_U238_0 = 0.70

    # Absorption cross-sections (in barns)
    sigma_a_U235 = 591.0
    sigma_a_U238 = 2.42
    sigma_a_Pu239 = 973.0

    # Fission cross-sections (in barns)
    sigma_f_U235 = 505.0
    sigma_f_Pu239 = 698.0

    # Total uranium burnup fraction
    burnup = 0.35

    print("Step-by-step calculation:")
    print("--------------------------")
    print("This script calculates the power fraction from Pu-239 after 35% uranium burnup.")
    print("Assumptions: Instantaneous decay of U-239/Np-239, one-group constant cross-sections.\n")

    # --- Step 1: Determine final uranium concentrations ---
    # Total uranium remaining is 1 - burnup
    N_U_total_final = 1.0 - burnup
    # Due to the much larger cross-section of U-235, it depletes much faster.
    # We can approximate that the remaining uranium is almost all U-238.
    N_U238_final = N_U_total_final
    print(f"1. Initial Conditions: N_U235_0 = {N_U235_0}, N_U238_0 = {N_U238_0}")
    print(f"2. After 35% burnup, total remaining Uranium is {N_U_total_final:.2f}.")
    print(f"   Approximating that U-235 is depleted, N_U238_final ≈ {N_U238_final:.4f}\n")

    # --- Step 2: Calculate the required neutron fluence (tau) ---
    # N_final = N_initial * exp(-sigma_a * tau) => tau = -ln(N_final / N_initial) / sigma_a
    tau_final = -math.log(N_U238_final / N_U238_0) / sigma_a_U238
    print(f"3. Calculate neutron fluence (τ) to reach this U-238 level:")
    print(f"   τ = -ln({N_U238_final:.4f} / {N_U238_0}) / {sigma_a_U238} = {tau_final:.6f} [neutrons/barn]\n")

    # --- Step 3: Calculate precise final nuclide concentrations using the fluence ---
    # Calculate the tiny amount of U-235 remaining
    N_U235_final = N_U235_0 * math.exp(-sigma_a_U235 * tau_final)

    # Calculate the amount of Pu-239 using the Bateman equation solution for the U-238 -> Pu-239 chain
    # N_Pu(t) = [N_U238_0 * σ_a_U238 / (σ_a_Pu239 - σ_a_U238)] * [exp(-σ_a_U238*τ) - exp(-σ_a_Pu239*τ)]
    prefactor = (N_U238_0 * sigma_a_U238) / (sigma_a_Pu239 - sigma_a_U238)
    exp_term_U238 = math.exp(-sigma_a_U238 * tau_final)
    exp_term_Pu239 = math.exp(-sigma_a_Pu239 * tau_final)
    N_Pu239_final = prefactor * (exp_term_U238 - exp_term_Pu239)

    print("4. Using this fluence, calculate the final concentrations of U-235 and Pu-239:")
    print(f"   N_U235_final = {N_U235_0} * exp(-{sigma_a_U235} * {tau_final:.6f}) = {N_U235_final:.4e}")
    print(f"   N_Pu239_final = {N_Pu239_final:.6f}\n")

    # --- Step 4: Calculate the power fraction from Pu-239 ---
    # Power is proportional to Fission Rate (N * sigma_f)
    power_from_U235 = N_U235_final * sigma_f_U235
    power_from_Pu239 = N_Pu239_final * sigma_f_Pu239
    total_power = power_from_U235 + power_from_Pu239
    fraction_from_Pu239 = power_from_Pu239 / total_power

    print("5. Calculate the power fraction from Plutonium-239.")
    print("   Power is proportional to (Number of Atoms * Fission Cross-Section).")
    print("   Fraction = (Power from Pu-239) / (Total Power)\n")
    
    print("Final Equation:")
    print("-------------")
    print(f"Fraction = (N_Pu239_final * σ_f_Pu239) / ( (N_U235_final * σ_f_U235) + (N_Pu239_final * σ_f_Pu239) )")
    print(f"Fraction = ({N_Pu239_final:.6f} * {sigma_f_Pu239}) / ( ({N_U235_final:.4e} * {sigma_f_U235}) + ({N_Pu239_final:.6f} * {sigma_f_Pu239}) )")
    print(f"Fraction = {power_from_Pu239:.6f} / ( {power_from_U235:.6e} + {power_from_Pu239:.6f} )")
    print(f"Fraction = {power_from_Pu239:.6f} / {total_power:.6f}\n")

    print("Final Answer:")
    print("-------------")
    print(f"The fraction of the power being produced from plutonium-239 is: {fraction_from_Pu239:.6f}")
    
    # Final answer in the required format
    print(f"\n<<<{fraction_from_Pu239:.6f}>>>")

if __name__ == '__main__':
    solve_power_fraction()