import math

def check_exoplanet_density_answer():
    """
    Checks the correctness of the answer to the exoplanet density question.

    The function models the density of rocky planets based on their mass relative to Earth.
    For planets of the same composition, density (ρ) scales with mass (M) due to
    gravitational compression. A common approximation is R ∝ M^k, where R is the radius
    and k ≈ 0.27 for rocky planets.
    Since density ρ = M/V and V ∝ R³, it follows that ρ ∝ M / (M^k)³ = M^(1-3k).

    This model is used to calculate and compare the densities of the planets described.
    """
    # The final answer from the LLM analysis to be checked.
    provided_answer = "C"

    # --- Model Parameters ---
    # Baseline density of Earth in g/cm^3
    earth_density = 5.51
    # Empirical exponent for the mass-radius relationship of rocky planets (R ∝ M^k)
    k = 0.27
    # Derived exponent for the mass-density relationship (ρ ∝ M^(1-3k))
    density_exponent = 1 - 3 * k

    # --- Calculate Densities for Each Option ---
    densities = {}

    # Option a: An Earth-mass and Earth-radius planet (our baseline).
    densities['a'] = earth_density

    # Option b: A planet with 2 Earth masses and a density of approximately 5.5 g/cm^3.
    # The density is explicitly given.
    densities['b'] = 5.5

    # Option c: A planet with the same composition as Earth but 5 times more massive.
    mass_ratio_c = 5.0
    densities['c'] = earth_density * (mass_ratio_c ** density_exponent)

    # Option d: A planet with the same composition as Earth but half the mass of Earth.
    mass_ratio_d = 0.5
    densities['d'] = earth_density * (mass_ratio_d ** density_exponent)

    # --- Determine the Correct Answer ---
    # Find the option with the highest calculated density.
    highest_density_option_id = max(densities, key=densities.get)

    # The question maps the descriptions (a, b, c, d) to letter choices (A, B, C, D).
    correct_letter_choice = highest_density_option_id.upper()

    # --- Verify and Return Result ---
    if provided_answer == correct_letter_choice:
        return "Correct"
    else:
        error_message = (
            f"Incorrect. The provided answer is '{provided_answer}', but the correct answer should be '{correct_letter_choice}'.\n\n"
            "Reasoning:\n"
            "The fundamental principle is that for planets of the same composition, greater mass leads to stronger gravitational compression and thus higher density.\n\n"
            "A quantitative check confirms this:\n"
            f"1. Density of planet (a) (Earth baseline): {densities['a']:.2f} g/cm³\n"
            f"2. Density of planet (b) (given): {densities['b']:.2f} g/cm³\n"
            f"3. Density of planet (c) (5x Earth mass): {densities['c']:.2f} g/cm³ (calculated based on compression)\n"
            f"4. Density of planet (d) (0.5x Earth mass): {densities['d']:.2f} g/cm³ (calculated based on compression)\n\n"
            f"Comparing these values, planet (c) has the highest density. This corresponds to letter choice '{correct_letter_choice}'."
        )
        return error_message

# Execute the check and print the result.
result = check_exoplanet_density_answer()
print(result)