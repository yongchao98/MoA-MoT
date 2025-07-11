def calculate_neff_scenario():
    """
    Demonstrates the effect of a new particle decaying into neutrinos on N_eff.

    N_eff parameterizes the radiation energy density relative to photons.
    rho_radiation = rho_photon * (1 + K * N_eff), where K is a constant.
    This implies rho_neutrino = rho_photon * K * N_eff.
    Therefore, N_eff = rho_neutrino / (rho_photon * K).
    """

    # --- 1. Define Standard Model parameters ---
    # N_eff in the Standard Model (slightly > 3 due to non-instantaneous decoupling)
    N_eff_SM = 3.044
    
    # The constant K combines statistical factors (7/8) and temperature differences.
    # K = (7/8) * (T_nu / T_gamma)^4. After e+e- annihilation, T_nu/T_gamma = (4/11)^(1/3).
    K = (7/8) * (4/11)**(4/3)

    # Let's assume a normalized photon energy density for this example.
    rho_photon = 100.0

    # Calculate the neutrino energy density in the Standard Model.
    rho_neutrino_SM = rho_photon * K * N_eff_SM

    print("--- Standard Cosmological Model ---")
    print(f"Assumed Photon Energy Density (rho_photon): {rho_photon:.2f}")
    print(f"Standard Model N_eff (N_eff_SM): {N_eff_SM}")
    print(f"Calculated SM Neutrino Energy Density (rho_neutrino_SM): {rho_neutrino_SM:.2f}\n")

    # --- 2. Introduce the New Physics Particle Decay ---
    # This particle decays, injecting extra energy into the neutrino population.
    # Let's assume this adds an arbitrary amount of energy density.
    delta_rho_from_decay = 15.0 # Energy density from new particle decay

    # The new total neutrino energy density is the sum of the standard one and the new injection.
    rho_neutrino_new = rho_neutrino_SM + delta_rho_from_decay
    
    print("--- New Physics Scenario ---")
    print(f"Energy injected into neutrinos from decay (delta_rho): {delta_rho_from_decay:.2f}")
    print(f"New total neutrino energy density (rho_neutrino_new): {rho_neutrino_SM:.2f} + {delta_rho_from_decay:.2f} = {rho_neutrino_new:.2f}\n")
    
    # --- 3. Calculate the new N_eff ---
    # The new N_eff is calculated from the new neutrino energy density.
    # N_eff_new = rho_neutrino_new / (rho_photon * K)
    N_eff_new = rho_neutrino_new / (rho_photon * K)

    print("--- Final Calculation ---")
    print("The new N_eff is calculated from the new total neutrino energy density.")
    print("Equation: N_eff_new = rho_neutrino_new / (rho_photon * K)")
    print(f"Calculation: {N_eff_new:.3f} = {rho_neutrino_new:.2f} / ({rho_photon:.2f} * {K:.4f})")
    
    # --- 4. Conclusion ---
    print("\n--- Conclusion ---")
    print(f"The new value N_eff_new ({N_eff_new:.3f}) is greater than the Standard Model value N_eff_SM ({N_eff_SM:.3f}).")
    print("Therefore, the decay of a new particle into neutrinos would increase N_eff.")

if __name__ == '__main__':
    calculate_neff_scenario()