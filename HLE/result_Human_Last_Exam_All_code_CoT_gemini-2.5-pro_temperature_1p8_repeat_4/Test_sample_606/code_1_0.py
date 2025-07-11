import sys

def solve_neff_change():
    """
    Analyzes the effect on N_eff from a new particle decaying into neutrinos.

    The function explains the physics and demonstrates that N_eff increases
    by showing a sample calculation.
    """

    # In the Standard Model (SM), N_eff accounts for 3 neutrino species, with
    # small corrections bringing the value to ~3.044.
    N_eff_SM = 3.044

    # Let's represent the SM neutrino energy density by a relative value of 1.
    # N_eff is directly proportional to this energy density.
    rho_nu_SM = 1.0

    # The new particle decays, adding a non-negligible amount of energy to the
    # neutrino sector. For this example, let's assume it increases the neutrino
    # energy density by a relative amount of 0.5. This value must be positive.
    rho_decay = 0.5

    # The new total energy density in the neutrino sector is the sum of the
    # standard contribution and the new decay energy.
    rho_nu_new = rho_nu_SM + rho_decay

    # The new N_eff will be scaled proportionally to the increase in energy density.
    N_eff_new = N_eff_SM * (rho_nu_new / rho_nu_SM)

    # The change, or delta, in N_eff is the difference.
    delta_N_eff = N_eff_new - N_eff_SM

    # --- Output the Explanation and Result ---
    print("Analysis of N_eff:")
    print("----------------")
    print(f"The Standard Model value, N_eff_SM, is ~{N_eff_SM:.3f}.")
    print("This value is proportional to the energy density of the standard neutrinos.")
    print("\nThe new particle decays exclusively into neutrinos, injecting extra energy into them.")
    print("This increases the total energy density of the neutrino sector while leaving the photon sector unaffected.")
    print("\nTherefore, N_eff must increase. Let's demonstrate with a calculation:")

    # Print the final equation as requested.
    print(f"\nNew N_eff = N_eff_SM + Delta_N_eff")
    print(f"{N_eff_new:.3f} = {N_eff_SM:.3f} + {delta_N_eff:.3f}")
    
    # A final clear statement of the outcome.
    print(f"\nSince {N_eff_new:.3f} > {N_eff_SM:.3f}, Neff increases.")


if __name__ == "__main__":
    solve_neff_change()
