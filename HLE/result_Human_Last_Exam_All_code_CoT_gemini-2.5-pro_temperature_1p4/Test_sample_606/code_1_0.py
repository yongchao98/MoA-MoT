def explain_neff_change():
    """
    Explains whether Neff would increase or decrease due to a hypothetical
    particle decaying into neutrinos in the early universe.
    """
    # 1. Define Neff and the Standard Model value.
    # Neff parameterizes the total relativistic energy density (radiation)
    # relative to the photon energy density.
    N_eff_SM = 3.044
    print(f"In the Standard Model of cosmology, the effective number of neutrino species is N_eff_SM = {N_eff_SM}.")
    print("This value corresponds to the energy density of the three known neutrino species.")
    print("-" * 50)

    # 2. Relationship between energy densities.
    # The total neutrino energy density (rho_nu) is directly proportional to Neff.
    # Let's represent the Standard Model energy densities symbolically.
    # rho_nu_SM is the standard neutrino energy density.
    # delta_rho_nu is the extra energy density from the new particle's decay.
    print("The energy density of neutrinos, rho_nu, is directly proportional to N_eff.")
    print("So, we can write: rho_nu_SM ∝ N_eff_SM")
    print("-" * 50)

    # 3. Introduce the new particle and its decay.
    # The new particle has some energy density, rho_X.
    # It decays entirely into neutrinos, injecting energy into them.
    print("A new, out-of-equilibrium particle exists and decays exclusively into neutrinos.")
    print("This decay injects an additional energy density, Δρ_nu, into the neutrino bath.")
    print("Since the particle has a non-negligible abundance, Δρ_nu > 0.")
    print("-" * 50)

    # 4. Calculate the new total neutrino energy density.
    # The new energy density is the sum of the old one and the injected energy.
    print("The new total neutrino energy density, rho_nu_new, is the sum of the standard amount and the new contribution:")
    print("rho_nu_new = rho_nu_SM + Δρ_nu")
    print("-" * 50)

    # 5. Determine the effect on Neff.
    # Since rho_nu_new > rho_nu_SM, the new Neff must be larger.
    print("Because rho_nu_new is greater than rho_nu_SM, the new value for N_eff, let's call it N_eff_new, must be greater than the standard value.")
    print("The final relationship can be expressed as an equation:")

    # Using placeholders for the symbolic equation
    delta_N_eff = "ΔN_eff" # This represents the positive increase
    print(f"\nFinal Equation: N_eff_new = N_eff_SM + {delta_N_eff}")
    print(f"Which means: N_eff_new = {N_eff_SM} + {delta_N_eff}, where {delta_N_eff} > 0.")

    print("\nTherefore, Neff would increase compared to the standard cosmological model.")

if __name__ == "__main__":
    explain_neff_change()