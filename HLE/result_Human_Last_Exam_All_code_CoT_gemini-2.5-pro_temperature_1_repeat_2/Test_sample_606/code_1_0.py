def calculate_neff_change():
    """
    Illustrates the effect on N_eff from a new particle decaying into neutrinos.

    This function models the change in the effective number of neutrino species (N_eff)
    when a hypothetical particle injects energy into the neutrino background.
    """
    print("--- Step 1: Define the Standard Model (SM) baseline ---")

    # In the SM, there are 3 neutrino species.
    # We can normalize the energy density of a single neutrino species to 1 for simplicity.
    rho_one_neutrino_species = 1.0
    print(f"Energy density of one standard model neutrino species (normalized): {rho_one_neutrino_species}")

    # The SM N_eff is ~3.044. We'll use 3 for this simplified model.
    N_eff_SM = 3.0
    # The total neutrino energy density in the SM is N_eff * rho_one_neutrino_species.
    rho_neutrino_SM = N_eff_SM * rho_one_neutrino_species
    print(f"Standard Model N_eff: {N_eff_SM}")
    print(f"Total neutrino energy density in SM: {N_eff_SM} * {rho_one_neutrino_species} = {rho_neutrino_SM}")

    print("\n--- Step 2: Introduce the new particle's decay ---")
    # A new particle with non-negligible abundance decays solely into neutrinos.
    # This injects a positive amount of energy density into the neutrino sector.
    # Let's assign an arbitrary positive value to this energy injection.
    injected_energy_density = 0.5
    print(f"Energy density injected into neutrinos from new particle decay: {injected_energy_density}")

    print("\n--- Step 3: Calculate the new total neutrino energy density ---")
    # The new total neutrino energy density is the SM value plus the injected amount.
    rho_neutrino_new = rho_neutrino_SM + injected_energy_density
    print(f"New total neutrino energy density = (SM density) + (injected density)")
    print(f"Equation: {rho_neutrino_new} = {rho_neutrino_SM} + {injected_energy_density}")


    print("\n--- Step 4: Calculate the new N_eff ---")
    # N_eff is defined as the total neutrino energy density divided by the
    # energy density of a single standard neutrino species.
    N_eff_new = rho_neutrino_new / rho_one_neutrino_species
    print(f"The new N_eff is calculated as (New total neutrino density) / (single species density)")
    print(f"Equation: {N_eff_new} = {rho_neutrino_new} / {rho_one_neutrino_species}")

    print("\n--- Step 5: Conclusion ---")
    if N_eff_new > N_eff_SM:
        result = "increase"
    elif N_eff_new < N_eff_SM:
        result = "decrease"
    else:
        result = "stay the same"

    print(f"The new value N_eff = {N_eff_new} is greater than the standard model value of {N_eff_SM}.")
    print(f"Therefore, N_eff would {result}.")


if __name__ == "__main__":
    calculate_neff_change()
