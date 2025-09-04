import collections

def check_exoplanet_density_answer():
    """
    Checks the correctness of the answer to the exoplanet density question.

    The function verifies the answer by:
    1. Establishing the physical principles for planetary density.
    2. Determining the correct planet option (a, b, c, or d).
    3. Applying the question's specific A/B/C/D mapping.
    4. Comparing the result with the provided answer.
    """
    # The final answer provided by the LLM being checked.
    llm_final_answer = "B"

    # 1. Define the planet options and the physical principles.
    # Principle: For planets of the same composition, density increases with mass
    # due to gravitational compression.
    # Let's represent the densities. We don't need exact values, just their order.
    # We can use a placeholder value for Earth's density.
    earth_density = 5.51  # g/cm^3

    planet_densities = {
        'a': earth_density,  # Baseline: Earth-mass and Earth-radius
        'b': 5.5,            # Given explicitly
        # 'c' is more massive with the same composition, so it's denser than Earth.
        'c': earth_density + 1, # Using a value > earth_density to represent this principle
        # 'd' is less massive with the same composition, so it's less dense than Earth.
        'd': earth_density - 1  # Using a value < earth_density to represent this principle
    }

    # 2. Identify the planet with the highest density based on our model.
    densest_planet = max(planet_densities, key=planet_densities.get)

    # 3. Define the mapping from the final answer choices (A, B, C, D) to the planet options.
    answer_mapping = {
        'A': 'd',
        'B': 'c',
        'C': 'a',
        'D': 'b'
    }

    # Find the correct final letter choice based on the densest planet.
    correct_final_answer = None
    for letter, planet in answer_mapping.items():
        if planet == densest_planet:
            correct_final_answer = letter
            break

    # 4. Check if the LLM's answer is correct.
    if llm_final_answer == correct_final_answer:
        return "Correct"
    else:
        reason = (
            f"Incorrect. The reasoning shows that planet '{densest_planet}' has the highest density. "
            f"According to the question's mapping (A->d, B->c, C->a, D->b), the correct final answer should be '{correct_final_answer}'. "
            f"The provided answer was '{llm_final_answer}', which corresponds to planet '{answer_mapping.get(llm_final_answer)}'."
        )
        return reason

# Execute the check and print the result.
result = check_exoplanet_density_answer()
print(result)