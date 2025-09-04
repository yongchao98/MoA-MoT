import math

def check_exoplanet_density_answer():
    """
    Checks the correctness of the answer to the exoplanet density question.

    The function models the density of four exoplanets based on the principles
    of gravitational compression and compares them to find the densest one.
    It then verifies if the provided answer corresponds to this result.
    """

    # --- Step 1: Define planet properties and calculate their densities ---

    # Baseline density of Earth in g/cm^3
    earth_density = 5.51

    # a) An Earth-mass and Earth-radius planet.
    # By definition, its density is Earth's density.
    density_a = earth_density

    # b) A planet with 2 Earth masses and a density of approximately 5.5 g/cm^3.
    # The density is explicitly given.
    density_b = 5.5

    # For planets (c) and (d), we use a simplified mass-density relationship
    # for rocky planets: ρ_planet = ρ_earth * (M_planet / M_earth)^α
    # A positive exponent α (derived from R ∝ M^k where k < 1/3) indicates
    # that density increases with mass due to gravitational compression.
    # A typical value for α is around 0.2.
    alpha = 0.2

    # c) A planet with the same composition as Earth but 5 times more massive.
    mass_ratio_c = 5.0
    density_c = earth_density * (mass_ratio_c ** alpha)

    # d) A planet with the same composition as Earth but half the mass of Earth.
    mass_ratio_d = 0.5
    density_d = earth_density * (mass_ratio_d ** alpha)

    # Store densities in a dictionary for easy comparison
    densities = {
        'a': density_a,
        'b': density_b,
        'c': density_c,
        'd': density_d
    }

    # --- Step 2: Identify the densest planet ---
    densest_planet_option = max(densities, key=densities.get)

    # --- Step 3: Verify the provided answer ---

    # The question maps the final answer letters (A, B, C, D) to the planet options (a, b, c, d)
    answer_mapping = {
        'A': 'a',
        'B': 'c',
        'C': 'd',
        'D': 'b'
    }

    # The final answer provided by the LLM
    provided_answer = 'B'

    # Check if the provided answer letter corresponds to the densest planet
    if answer_mapping.get(provided_answer) == densest_planet_option:
        return "Correct"
    else:
        # Find the correct answer letter for the densest planet
        correct_answer_letter = [key for key, val in answer_mapping.items() if val == densest_planet_option][0]
        
        reason = (
            f"The answer is incorrect. The analysis shows that planet '{densest_planet_option}' has the highest density.\n"
            f"Calculated densities (g/cm^3): a={densities['a']:.2f}, b={densities['b']:.2f}, c={densities['c']:.2f}, d={densities['d']:.2f}.\n"
            f"The densest planet is '{densest_planet_option}', which corresponds to answer choice '{correct_answer_letter}'.\n"
            f"The provided answer was '{provided_answer}', which corresponds to planet '{answer_mapping.get(provided_answer)}'."
        )
        return reason

# Execute the check and print the result
result = check_exoplanet_density_answer()
print(result)