import sys

# Set encoding to UTF-8 to support special characters like μ and ε₀
if sys.stdout.encoding != 'utf-8':
    sys.stdout.reconfigure(encoding='utf-8')

def explain_chemical_potential_limit():
    """
    Explains the fundamental limit on the chemical potential for bosons
    in Bose-Einstein condensation.
    """
    ground_state_energy = "ε₀"
    chemical_potential = "μ"
    boltzmann_const = "k_B"
    temperature = "T"
    energy_level_i = "ε_i"
    occupation_number_i = "<n_i>"

    print("1. The Bose-Einstein Distribution")
    print("In quantum statistics, the average number of bosons in a single-particle energy state 'i' is given by:")
    print(f"  {occupation_number_i} = 1 / (exp[({energy_level_i} - {chemical_potential}) / ({boltzmann_const}*{temperature})] - 1)")
    print("-" * 50)

    print("2. The Physical Constraint")
    print(f"For the occupation number {occupation_number_i} to be positive, the denominator must be positive.")
    print("This requires that for every energy level:")
    print(f"  {energy_level_i} - {chemical_potential} > 0")
    print(f"Which implies: {chemical_potential} < {energy_level_i}")
    print("\nThe most restrictive limit comes from the lowest possible energy state, the ground state, with energy " + ground_state_energy + ".")
    print("Therefore, the chemical potential must always be less than the ground state energy:")
    print(f"  {chemical_potential} < {ground_state_energy}")
    print("-" * 50)

    print("3. Behavior during Condensation")
    print("As a Bose gas is cooled towards the critical temperature (T_c), the chemical potential approaches the ground state energy from below.")
    print("For a system undergoing Bose-Einstein condensation (at T ≤ T_c), the chemical potential is effectively pinned to this limiting value to allow for macroscopic occupation of the ground state.")
    print(f"  {chemical_potential} = {ground_state_energy}  (for T ≤ T_c)")
    print("-" * 50)

    print("4. Conclusion")
    print("The chemical potential of a non-interacting Bose gas at absolute zero (T=0) is, by definition, the ground state energy, " + ground_state_energy + ".")
    print("Therefore, during Bose-Einstein condensation, the chemical potential must be equal to the chemical potential of a non-interacting Bose gas at zero temperature.")
    print("\nFinal statement of the limit:")
    print(f"'{chemical_potential} must be equal to the chemical potential of a non-interacting Bose gas at zero temperature.'")


if __name__ == "__main__":
    explain_chemical_potential_limit()