# Neff_calculation.py

def calculate_neff_change():
    """
    Calculates the change in N_eff due to a hypothetical particle decaying into neutrinos.
    """
    
    # 1. Theoretical background
    print("In the Standard Model of Cosmology, the total radiation energy density (rho_rad) is composed of photons (rho_gamma) and neutrinos (rho_nu).")
    print("This relationship is parameterized by N_eff, the effective number of neutrino species:")
    print("rho_rad = rho_gamma + rho_nu = rho_gamma * (1 + N_eff * (7/8) * (T_nu/T_gamma)^4)\n")

    # 2. Standard Model value
    N_eff_SM = 3.045
    print(f"The Standard Model predicts N_eff_SM = {N_eff_SM}.\n")

    # 3. Introduce the new physics
    print("We introduce a new massive particle that decays, adding energy (Δρ_nu) exclusively to the neutrino population.")
    print("The photon energy density (rho_gamma) is unaffected.")
    print("The new neutrino energy density is rho_nu_new = rho_nu_SM + Δρ_nu.\n")
    
    # Let's assume the additional energy from the decay, Δρ_nu, is a certain fraction
    # of the photon energy density, rho_gamma. This is a model-dependent parameter.
    # We choose a non-negligible value for demonstration.
    delta_rho_nu_over_rho_gamma = 0.2
    
    print(f"Let's assume this additional energy is {delta_rho_nu_over_rho_gamma:.2f} times the photon energy density (Δρ_nu / ρ_gamma = {delta_rho_nu_over_rho_gamma}).\n")
    
    # 4. Calculate the change in N_eff
    # From the definition, the change in N_eff is related to the additional energy density.
    # Δρ_nu = ΔN_eff * (7/8) * (T_nu/T_gamma)^4 * rho_gamma
    # ΔN_eff = (Δρ_nu / rho_gamma) / ( (7/8) * (T_nu/T_gamma)^4 )
    
    # In the SM, entropy conservation after neutrino decoupling gives T_nu/T_gamma = (4/11)^(1/3).
    T_ratio_4 = (4.0 / 11.0)**(4.0 / 3.0)
    
    # Denominator term
    factor = (7.0 / 8.0) * T_ratio_4
    
    # Calculate the change in N_eff
    delta_N_eff = delta_rho_nu_over_rho_gamma / factor
    
    # Calculate the new N_eff
    N_eff_new = N_eff_SM + delta_N_eff

    # 5. Print the step-by-step calculation of the final equation
    print("--- Calculation ---")
    print("The final value, N_eff_new, is calculated from the equation: N_eff_new = N_eff_SM + ΔN_eff")
    print(f"where ΔN_eff = (Δρ_nu / ρ_gamma) / [ (7/8) * (T_nu/T_gamma)^4 ]")
    print("\nPlugging in the numbers:")
    print(f"The ratio of temperatures term is (T_nu/T_gamma)^4 = (4/11)^(4/3) ≈ {T_ratio_4:.5f}")
    print(f"The denominator is (7/8) * {T_ratio_4:.5f} ≈ {factor:.5f}")
    print(f"So, ΔN_eff = {delta_rho_nu_over_rho_gamma} / {factor:.5f} ≈ {delta_N_eff:.5f}\n")
    
    print("--- Final Equation ---")
    print("N_eff_new = N_eff_SM + ΔN_eff")
    # As requested, printing each number in the final equation:
    print(f"{N_eff_new:.5f} = {N_eff_SM} + {delta_N_eff:.5f}")
    
    print("\n--- Conclusion ---")
    print(f"Because the particle decay adds energy to the neutrinos, ΔN_eff is positive.")
    print(f"Therefore, N_eff increases from {N_eff_SM} to {N_eff_new:.5f}.")

if __name__ == "__main__":
    calculate_neff_change()
