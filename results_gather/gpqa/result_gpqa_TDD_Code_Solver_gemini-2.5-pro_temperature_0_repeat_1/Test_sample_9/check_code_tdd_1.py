import math

def check_answer_correctness():
    """
    This function checks the correctness of the provided answer to the exoplanet density question.

    The question asks to identify the exoplanet with the highest density among four options:
    a) An Earth-mass and Earth-radius planet.
    b) A planet with 2 Earth masses and a density of approximately 5.5 g/cm^3.
    c) A planet with the same composition as Earth but 5 times more massive than Earth.
    d) A planet with the same composition as Earth but half the mass of Earth.

    The multiple-choice options are: A) b, B) d, C) a, D) c.
    The provided answer to check is 'D'.
    """
    provided_answer = 'D'

    # --- Step 1: Calculate the density for each scenario ---

    # Use Earth's average density as a baseline.
    # Earth's density is approximately 5.51 g/cm^3.
    earth_density = 5.51

    # Scenario a: An Earth-mass and Earth-radius planet.
    # By definition, its density is Earth's density.
    density_a = earth_density

    # Scenario b: A planet with a given density of approximately 5.5 g/cm^3.
    # The density is explicitly stated.
    density_b = 5.5

    # For scenarios c and d, we use the mass-density relationship for rocky planets.
    # Due to gravitational self-compression, a more massive planet of the same
    # composition will be denser. The relationship is approximately ρ ∝ M^(1/4).
    
    # Scenario c: A planet with the same composition as Earth but 5 times more massive.
    mass_ratio_c = 5.0
    density_c = earth_density * (mass_ratio_c ** 0.25)

    # Scenario d: A planet with the same composition as Earth but half the mass.
    mass_ratio_d = 0.5
    density_d = earth_density * (mass_ratio_d ** 0.25)

    # --- Step 2: Identify the scenario with the highest density ---
    
    densities = {
        'a': density_a,
        'b': density_b,
        'c': density_c,
        'd': density_d
    }

    # Find the key (scenario letter) corresponding to the maximum density value.
    highest_density_scenario = max(densities, key=densities.get)

    # --- Step 3: Map the scenario to the correct multiple-choice answer ---
    
    # The question maps scenarios to choices as follows:
    # A) b, B) d, C) a, D) c
    answer_mapping = {
        'a': 'C',
        'b': 'A',
        'c': 'D',
        'd': 'B'
    }
    
    calculated_correct_answer = answer_mapping[highest_density_scenario]

    # --- Step 4: Compare the calculated answer with the provided answer ---

    if calculated_correct_answer == provided_answer:
        return "Correct"
    else:
        reason = (
            f"The provided answer '{provided_answer}' is incorrect.\n"
            f"The calculated densities (in g/cm^3) are:\n"
            f"  - Scenario 'a': {densities['a']:.2f}\n"
            f"  - Scenario 'b': {densities['b']:.2f}\n"
            f"  - Scenario 'c': {densities['c']:.2f}\n"
            f"  - Scenario 'd': {densities['d']:.2f}\n"
            f"The highest density belongs to scenario '{highest_density_scenario}' "
            f"with a value of {densities[highest_density_scenario]:.2f} g/cm^3.\n"
            f"According to the answer mapping, scenario '{highest_density_scenario}' corresponds to option '{calculated_correct_answer}'.\n"
            f"Therefore, the correct answer is '{calculated_correct_answer}'."
        )
        return reason

# Execute the check
result = check_answer_correctness()
print(result)