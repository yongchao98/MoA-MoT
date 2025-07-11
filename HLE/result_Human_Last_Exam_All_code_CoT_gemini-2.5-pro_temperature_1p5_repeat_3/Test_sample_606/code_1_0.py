import math

def calculate_neff_change():
    """
    Illustrates how an energy injection into the neutrino background affects Neff.
    """
    
    # 1. Standard Model Parameters
    # Neff in the Standard Model is slightly greater than 3 due to non-instantaneous
    # decoupling and QED finite temperature corrections.
    N_eff_SM = 3.044
    print(f"Standard Model Neff (N_eff_SM): {N_eff_SM}")
    print("-" * 30)

    # 2. Relationship between Neff and Energy Densities
    # The energy density of neutrinos (rho_nu) is related to the energy density of
    # photons (rho_gamma) by the formula:
    # rho_nu = N_eff * (7/8) * (4/11)^(4/3) * rho_gamma
    # Let's define the constant factor in this relation.
    
    factor = (7.0 / 8.0) * (4.0 / 11.0)**(4.0 / 3.0)
    print(f"The constant factor relating Neff to energy density ratio is: {factor:.4f}")

    # For this example, let's assume a normalized photon energy density.
    rho_gamma = 100.0  # Arbitrary units
    print(f"Assuming a photon energy density (rho_gamma): {rho_gamma}")
    
    # Calculate the neutrino energy density in the Standard Model.
    rho_nu_SM = N_eff_SM * factor * rho_gamma
    print(f"This gives a standard neutrino energy density (rho_nu_SM) = {N_eff_SM} * {factor:.4f} * {rho_gamma} = {rho_nu_SM:.2f}")
    print("-" * 30)
    
    # 3. Effect of the New Particle Decay
    # The new particle decays and injects energy solely into neutrinos.
    # Let's assume this adds 10 units of energy density.
    energy_injection_from_decay = 10.0
    print(f"A new particle decays, injecting energy into neutrinos: +{energy_injection_from_decay}")

    # Calculate the new total neutrino energy density.
    rho_nu_new = rho_nu_SM + energy_injection_from_decay
    print(f"The new neutrino energy density (rho_nu_new) = {rho_nu_SM:.2f} + {energy_injection_from_decay} = {rho_nu_new:.2f}")
    
    # 4. Calculate the New Neff
    # We can now rearrange the formula to solve for the new Neff.
    # Neff_new = rho_nu_new / (factor * rho_gamma)
    N_eff_new = rho_nu_new / (factor * rho_gamma)
    
    print("\nCalculating the new value of Neff:")
    print(f"N_eff_new = rho_nu_new / (factor * rho_gamma)")
    print(f"N_eff_new = {rho_nu_new:.2f} / ({factor:.4f} * {rho_gamma})")
    print(f"N_eff_new = {rho_nu_new:.2f} / {factor * rho_gamma:.2f}")
    print(f"N_eff_new = {N_eff_new:.4f}")
    print("-" * 30)

    # 5. Conclusion
    print(f"Final Comparison:")
    print(f"Initial Neff: {N_eff_SM}")
    print(f"Final Neff:   {N_eff_new:.4f}")
    if N_eff_new > N_eff_SM:
        print("\nConclusion: Neff has increased.")
    else:
        print("\nConclusion: Neff has not increased.")

if __name__ == "__main__":
    calculate_neff_change()
