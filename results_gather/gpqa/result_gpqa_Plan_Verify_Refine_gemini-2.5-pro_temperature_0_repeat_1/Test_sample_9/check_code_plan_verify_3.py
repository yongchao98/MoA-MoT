import math

def check_correctness_of_exoplanet_density_answer():
    """
    Checks the correctness of the answer to the exoplanet density question.

    The function models the density of rocky planets to verify the provided answer.
    The key physical principle is gravitational compression: for planets of the same
    composition, higher mass leads to higher density.

    This is modeled by the mass-radius relationship R ∝ M^β (where β ≈ 0.27 for rocky planets).
    From this, it follows that density ρ ∝ M^(1 - 3β). Since 1 - 3β > 0, density
    increases with mass.
    """
    llm_provided_answer = 'C'

    # --- Define constants and model parameters ---
    # Earth's average density in g/cm^3
    density_earth = 5.51

    # Exponent in the mass-radius relationship R ∝ M^β for rocky planets
    beta = 0.27
    # Resulting exponent in the mass-density relationship ρ ∝ M^(1 - 3β)
    density_exponent = 1 - 3 * beta

    # --- Calculate/assign density for each option ---
    densities = {}

    # a) An Earth-mass and Earth-radius planet (i.e., Earth). This is our baseline.
    densities['a'] = density_earth

    # b) A planet with a given density of approximately 5.5 g/cm^3.
    densities['b'] = 5.5

    # c) A planet with the same composition as Earth but 5 times more massive.
    # We use the mass-density relationship: ρ_planet = ρ_earth * (M_planet / M_earth)^exponent
    mass_ratio_c = 5.0
    densities['c'] = density_earth * (mass_ratio_c ** density_exponent)

    # d) A planet with the same composition as Earth but half the mass of Earth.
    mass_ratio_d = 0.5
    densities['d'] = density_earth * (mass_ratio_d ** density_exponent)

    # --- Verification ---
    # Find the option corresponding to the maximum calculated density
    calculated_highest_density_option = max(densities, key=densities.get)

    # Compare the calculated result with the LLM's answer
    if calculated_highest_density_option.upper() == llm_provided_answer:
        return "Correct"
    else:
        reason = (
            f"The answer is incorrect. The reasoning should lead to a different conclusion.\n"
            f"The principle of gravitational compression states that for planets of the same composition, "
            f"a more massive planet will be denser. Let's quantify this:\n\n"
            f"Calculated densities:\n"
            f"  - Option a: {densities['a']:.2f} g/cm³ (Earth baseline)\n"
            f"  - Option b: {densities['b']:.2f} g/cm³ (Given value)\n"
            f"  - Option c (5x Earth mass): {densities['c']:.2f} g/cm³\n"
            f"  - Option d (0.5x Earth mass): {densities['d']:.2f} g/cm³\n\n"
            f"The highest density belongs to option '{calculated_highest_density_option}', which is more massive than Earth. "
            f"The provided answer was '{llm_provided_answer}'."
        )
        return reason

# Run the check
result = check_correctness_of_exoplanet_density_answer()
print(result)