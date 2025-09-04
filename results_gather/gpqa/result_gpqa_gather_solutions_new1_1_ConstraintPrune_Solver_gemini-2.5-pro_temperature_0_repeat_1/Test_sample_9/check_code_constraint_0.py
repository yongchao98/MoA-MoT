import math

def check_exoplanet_density_answer():
    """
    Checks the correctness of the answer to the exoplanet density question.

    The function models the densities of the four planets based on physical principles
    and given information, identifies the densest planet, and verifies if it matches
    the provided answer choice according to the question's mapping.
    """

    # --- Constants and Model Parameters ---
    # Earth's average density in g/cm^3
    RHO_EARTH = 5.51

    # Mass-radius exponent 'k' for rocky planets (where R ∝ M^k).
    # This value is based on planetary models and is typically around 0.27.
    K_EXPONENT = 0.27

    # --- Calculate/Define Densities for each planet option ---
    densities = {}

    # Option a: An Earth-mass and Earth-radius planet.
    # By definition, its density is Earth's density.
    densities['a'] = RHO_EARTH

    # Option b: A planet with 2 Earth masses and a density of approximately 5.5 g/cm^3.
    # The density is explicitly given in the problem statement.
    densities['b'] = 5.5

    # Option c: A planet with the same composition as Earth but 5 times more massive.
    # We use the derived mass-density relationship: ρ_planet = ρ_earth * (M_ratio)^(1 - 3*k)
    mass_ratio_c = 5.0
    densities['c'] = RHO_EARTH * (mass_ratio_c ** (1 - 3 * K_EXPONENT))

    # Option d: A planet with the same composition as Earth but half the mass of Earth.
    mass_ratio_d = 0.5
    densities['d'] = RHO_EARTH * (mass_ratio_d ** (1 - 3 * K_EXPONENT))

    # --- Identify the densest planet from our calculations ---
    densest_planet_label = max(densities, key=densities.get)

    # --- Verify against the provided answer and mapping ---
    # The question provides a specific mapping from answer choice to planet option.
    answer_mapping = {
        'A': 'a',
        'B': 'c',
        'C': 'd',
        'D': 'b'
    }

    # The final answer provided by the LLM to be checked.
    provided_answer_choice = 'B'

    # Find which planet corresponds to the provided answer choice.
    planet_from_answer = answer_mapping.get(provided_answer_choice)

    # Check if the planet from the answer is the one we calculated as the densest.
    if planet_from_answer == densest_planet_label:
        return "Correct"
    else:
        # If the answer is incorrect, provide a detailed reason.
        reason = (
            f"The provided answer is '{provided_answer_choice}', which corresponds to planet '{planet_from_answer}'.\n"
            f"However, the physical model shows that planet '{densest_planet_label}' has the highest density.\n\n"
            f"Calculated densities (in g/cm^3):\n"
            f"  - Planet a: {densities['a']:.2f} (baseline)\n"
            f"  - Planet b: {densities['b']:.2f} (given)\n"
            f"  - Planet c: {densities['c']:.2f} (calculated)\n"
            f"  - Planet d: {densities['d']:.2f} (calculated)\n\n"
            f"The densest planet is '{densest_planet_label}' with a density of {densities[densest_planet_label]:.2f} g/cm^3. "
            f"Therefore, the correct answer choice should be the one that maps to '{densest_planet_label}'."
        )
        return reason

# Execute the check and print the result.
result = check_exoplanet_density_answer()
print(result)