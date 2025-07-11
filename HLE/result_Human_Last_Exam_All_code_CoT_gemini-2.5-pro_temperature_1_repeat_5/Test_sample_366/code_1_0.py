import numpy as np
from scipy.optimize import fsolve

def solve_power_fraction():
    """
    Calculates the fraction of power from Pu-239 in a nuclear reactor
    after a specified burnup.
    """
    # Step 1: Define initial conditions and physical constants.
    # ---------------------------------------------------------
    # Initial relative number of atoms based on 30% enrichment.
    N_U235_0 = 0.30
    N_U238_0 = 0.70
    N_U_total_0 = N_U235_0 + N_U238_0

    # Microscopic cross-sections in barns (1 barn = 10^-24 cm^2).
    sigma_a_U235 = 591.0   # Absorption for U-235
    sigma_f_U235 = 505.0   # Fission for U-235
    sigma_a_U238 = 2.42    # Absorption for U-238
    sigma_a_Pu239 = 973.0  # Absorption for Pu-239
    sigma_f_Pu239 = 698.0  # Fission for Pu-239

    # Total uranium burnup fraction.
    burnup_fraction = 0.35

    print("This problem calculates the fraction of power produced by plutonium-239 after 35% of the initial uranium has been consumed.")
    print("\n--- Step 1: Initial Conditions and Constants ---")
    print(f"Initial U-235 fraction (N_U235_0): {N_U235_0}")
    print(f"Initial U-238 fraction (N_U238_0): {N_U238_0}")
    print(f"σ_a for U-235: {sigma_a_U235} barns, σ_f for U-235: {sigma_f_U235} barns")
    print(f"σ_a for U-238: {sigma_a_U238} barns")
    print(f"σ_a for Pu-239: {sigma_a_Pu239} barns, σ_f for Pu-239: {sigma_f_Pu239} barns")
    print(f"Uranium burnup fraction: {burnup_fraction}")

    # Step 2: Calculate the neutron fluence (Φ) for the given burnup.
    # ----------------------------------------------------------------
    # We solve the equation: N_U_total(Φ) = (1 - burnup_fraction) * N_U_total_0
    # N_U235_0*exp(-σ_a_U235*Φ) + N_U238_0*exp(-σ_a_U238*Φ) = (1 - 0.35) * 1.0
    def find_fluence_eq(fluence):
        """Equation whose root is the required neutron fluence."""
        n_u235_t = N_U235_0 * np.exp(-sigma_a_U235 * fluence)
        n_u238_t = N_U238_0 * np.exp(-sigma_a_U238 * fluence)
        target_n_u_total = (1 - burnup_fraction) * N_U_total_0
        return n_u235_t + n_u238_t - target_n_u_total

    # Solve for fluence using a numerical solver.
    initial_guess_fluence = 0.03 # An educated guess to help the solver.
    fluence = fsolve(find_fluence_eq, initial_guess_fluence)[0]

    print("\n--- Step 2: Calculate Neutron Fluence (Φ) ---")
    print("We solve the equation: N_U235_0*exp(-σ_a_U235*Φ) + N_U238_0*exp(-σ_a_U238*Φ) = 1 - burnup_fraction")
    print(f"The calculated neutron fluence (Φ) is: {fluence:.6f} [1/barns]")

    # Step 3: Calculate the final concentrations of U-235 and Pu-239.
    # -----------------------------------------------------------------
    # Equation for U-235 depletion.
    N_U235_final = N_U235_0 * np.exp(-sigma_a_U235 * fluence)

    # Equation for Pu-239 production and depletion (solution to Bateman equation).
    term_pu_prefactor = (N_U238_0 * sigma_a_U238) / (sigma_a_Pu239 - sigma_a_U238)
    term_pu_exp_diff = np.exp(-sigma_a_U238 * fluence) - np.exp(-sigma_a_Pu239 * fluence)
    N_Pu239_final = term_pu_prefactor * term_pu_exp_diff

    print("\n--- Step 3: Calculate Final Nuclide Concentrations ---")
    print("Using the fluence, we find the final relative number of atoms (N_final) for U-235 and Pu-239.")
    print(f"Final relative concentration of U-235 (N_U235_final): {N_U235_final:.3e}")
    print(f"Final relative concentration of Pu-239 (N_Pu239_final): {N_Pu239_final:.6f}")

    # Step 4: Calculate the power fraction from Pu-239.
    # ----------------------------------------------------
    # Power is proportional to the macroscopic fission rate (R_f = N * σ_f * flux).
    # The flux term cancels in the ratio.
    R_f_U235 = N_U235_final * sigma_f_U235
    R_f_Pu239 = N_Pu239_final * sigma_f_Pu239
    
    # Total fission rate is the sum of individual rates.
    R_f_total = R_f_U235 + R_f_Pu239

    # The fraction of power from Pu-239.
    power_fraction_pu = R_f_Pu239 / R_f_total

    print("\n--- Step 4: Calculate Power Fraction from Plutonium-239 ---")
    print("Power is proportional to the fission rate (R_f = N * σ_f). The fraction from Pu-239 is R_f_Pu239 / (R_f_U235 + R_f_Pu239).")
    print("\nFinal Equation:")
    print(f"Fraction = ({N_Pu239_final:.6f} * {sigma_f_Pu239}) / (({N_U235_final:.3e} * {sigma_f_U235}) + ({N_Pu239_final:.6f} * {sigma_f_Pu239}))")
    print(f"          = {R_f_Pu239:.6f} / ({R_f_U235:.3e} + {R_f_Pu239:.6f})")
    
    print(f"\nResult:")
    print(f"The fraction of power produced from plutonium-239 is: {power_fraction_pu:.6f}")
    
    print(f"\n<<< {power_fraction_pu:.6f} >>>")

if __name__ == '__main__':
    solve_power_fraction()