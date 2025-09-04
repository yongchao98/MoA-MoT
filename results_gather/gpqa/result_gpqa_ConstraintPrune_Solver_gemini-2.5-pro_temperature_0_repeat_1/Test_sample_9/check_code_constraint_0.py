import math

def check_exoplanet_density():
    """
    This function checks which of the described exoplanets has the highest density.
    It uses the physical principle that for a given composition, density increases
    with mass due to gravitational self-compression.
    """

    # --- Define Baseline and Planet Properties ---

    # Use Earth's density as a baseline.
    # The exact value isn't critical, as we are looking for the *highest* density.
    earth_density = 5.51  # in g/cm^3

    # Option a: An Earth-mass and Earth-radius planet.
    # Its density is by definition Earth's density.
    density_a = earth_density

    # Option b: A planet with 2 Earth masses and a density of approximately 5.5 g/cm^3.
    # The density is explicitly given.
    density_b = 5.5

    # Options c and d: Same composition as Earth, different masses.
    # For rocky planets, a common mass-radius relationship is R ∝ M^α, where α is ~0.25-0.3.
    # Density (ρ) = Mass / Volume ∝ M / R^3 ∝ M / (M^α)^3 ∝ M / M^(3α) ∝ M^(1-3α).
    # Since α < 1/3, the exponent (1-3α) is positive, meaning density increases with mass.
    # We'll use α = 0.27 as a reasonable exponent for rocky planets.
    # So, ρ_planet = ρ_earth * (M_planet / M_earth)^(1 - 3*0.27)
    exponent = 1 - 3 * 0.27 # This is approx 0.19

    # Option c: 5 times more massive than Earth.
    mass_ratio_c = 5.0
    density_c = earth_density * (mass_ratio_c ** exponent)

    # Option d: Half the mass of Earth.
    mass_ratio_d = 0.5
    density_d = earth_density * (mass_ratio_d ** exponent)

    # --- Verification ---

    # Store the calculated densities in a dictionary for easy comparison.
    densities = {
        "a": density_a,
        "b": density_b,
        "c": density_c,
        "d": density_d
    }

    # Find the option with the highest calculated density.
    calculated_winner = max(densities, key=densities.get)

    # The provided answer is 'A', which corresponds to option 'c' in the question.
    # Let's define the mapping from the answer choices (A, B, C, D) to the planet descriptions (a, b, c, d).
    answer_mapping = {'A': 'c', 'B': 'b', 'C': 'd', 'D': 'a'}
    provided_answer_choice = 'A'
    llm_choice = answer_mapping.get(provided_answer_choice)

    # Check if the LLM's choice matches the calculated winner.
    if llm_choice == calculated_winner:
        return "Correct"
    else:
        reason = f"The answer is incorrect.\n"
        reason += f"The provided answer is '{provided_answer_choice}', which corresponds to planet '{llm_choice}'.\n"
        reason += "However, based on physical principles, the planet with the highest density is '{calculated_winner}'.\n"
        reason += "Calculated densities:\n"
        reason += f"  - Planet a: {densities['a']:.2f} g/cm^3 (Earth baseline)\n"
        reason += f"  - Planet b: {densities['b']:.2f} g/cm^3 (as stated in the option)\n"
        reason += f"  - Planet c: {densities['c']:.2f} g/cm^3 (more massive, so more compressed and denser)\n"
        reason += f"  - Planet d: {densities['d']:.2f} g/cm^3 (less massive, so less compressed and less dense)\n"
        reason += f"The highest density belongs to planet '{calculated_winner}', not '{llm_choice}'."
        return reason

# Execute the check and print the result.
result = check_exoplanet_density()
print(result)