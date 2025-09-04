import math

def check_exoplanet_density_answer():
    """
    Checks the correctness of the answer for the exoplanet density question.

    The core principle is that for planets of the same composition, density increases
    with mass due to gravitational compression. A common approximation is that
    density (ρ) scales with mass (M) as ρ ∝ M^(1/4).

    The function calculates the density for each scenario, identifies the highest,
    and compares the result with the provided answer.
    """
    # The answer provided by the LLM's code.
    llm_answer = 'D'

    # --- Calculation ---
    # Constants and initial values
    earth_density = 5.51  # g/cm^3

    # a) An Earth-mass and Earth-radius planet.
    density_a = earth_density

    # b) A planet with a density of approximately 5.5 g/cm^3.
    density_b = 5.5

    # c) A planet with the same composition as Earth but 5 times more massive.
    # Using the relationship: ρ_new = ρ_earth * (M_new / M_earth)^(1/4)
    mass_ratio_c = 5.0
    density_c = earth_density * (mass_ratio_c ** 0.25)

    # d) A planet with the same composition as Earth but half the mass of Earth.
    mass_ratio_d = 0.5
    density_d = earth_density * (mass_ratio_d ** 0.25)

    # Store densities in a dictionary to find the maximum
    densities = {
        'a': density_a,
        'b': density_b,
        'c': density_c,
        'd': density_d
    }

    # Find the scenario with the highest density
    highest_density_scenario = max(densities, key=densities.get)

    # Map the scenario ('a', 'b', 'c', 'd') to the multiple-choice option ('A', 'B', 'C', 'D')
    # From the question: A) b, B) d, C) a, D) c
    answer_mapping = {'a': 'C', 'b': 'A', 'c': 'D', 'd': 'B'}
    correct_answer_option = answer_mapping[highest_density_scenario]

    # --- Verification ---
    if correct_answer_option == llm_answer:
        return "Correct"
    else:
        reason = f"The provided answer '{llm_answer}' is incorrect.\n"
        reason += "Calculated densities (in g/cm^3):\n"
        reason += f"  a) Earth-like: {densities['a']:.2f}\n"
        reason += f"  b) Given density: {densities['b']:.2f}\n"
        reason += f"  c) 5x Earth mass: {densities['c']:.2f}\n"
        reason += f"  d) 0.5x Earth mass: {densities['d']:.2f}\n"
        reason += f"The scenario with the highest density is '{highest_density_scenario}' ({densities[highest_density_scenario]:.2f} g/cm^3).\n"
        reason += f"This corresponds to multiple-choice option '{correct_answer_option}', not '{llm_answer}'."
        return reason

# Run the check
result = check_exoplanet_density_answer()
print(result)