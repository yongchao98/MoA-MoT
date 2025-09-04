import math

def check_exoplanet_density_answer():
    """
    Checks the correctness of the answer to the exoplanet density question.

    The core principle is that for rocky planets of the same composition,
    higher mass leads to greater gravitational compression and thus higher density.
    This can be modeled with a mass-radius relationship R ∝ M^k, where k < 1/3.
    This implies that density ρ ∝ M^(1-3k), so density increases with mass.
    """
    # The provided answer from the LLM is 'C'.
    llm_answer = 'C'

    # --- Define Constants and Model ---
    # We can work with ratios relative to Earth.
    # Earth's density in g/cm^3.
    rho_earth = 5.51

    # A reasonable exponent for the mass-radius relationship (R ∝ M^k) for rocky planets.
    # The exact value isn't critical, as long as k < 1/3, density will increase with mass.
    # A value of k=0.27 is a common approximation.
    mass_radius_exponent = 0.27

    # --- Calculate Density for Each Option ---

    # Option a: An Earth-mass and Earth-radius planet (our baseline).
    density_a = rho_earth

    # Option b: A planet with a given density of ~5.5 g/cm^3.
    density_b = 5.5

    # Option c: A planet with the same composition as Earth but 5x the mass.
    # Density (ρ) is proportional to Mass^(1 - 3*k).
    mass_ratio_c = 5.0
    density_c = rho_earth * (mass_ratio_c ** (1 - 3 * mass_radius_exponent))

    # Option d: A planet with the same composition as Earth but 0.5x the mass.
    mass_ratio_d = 0.5
    density_d = rho_earth * (mass_ratio_d ** (1 - 3 * mass_radius_exponent))

    # --- Compare Densities and Verify the Answer ---
    densities = {
        'A': density_a,
        'B': density_b,
        'C': density_c,
        'D': density_d
    }

    # Find the option with the highest calculated density.
    highest_density_option = max(densities, key=densities.get)

    # Check if the calculated correct option matches the LLM's answer.
    if highest_density_option == llm_answer:
        return "Correct"
    else:
        reason = (f"The answer '{llm_answer}' is incorrect.\n"
                  f"The analysis of densities shows that option '{highest_density_option}' has the highest density.\n"
                  f"Reasoning:\n"
                  f" - a) Earth's density is ~{densities['A']:.2f} g/cm³.\n"
                  f" - b) Density is given as ~{densities['B']:.2f} g/cm³.\n"
                  f" - c) A 5-Earth-mass planet with Earth's composition would be significantly compressed by its own gravity, "
                  f"leading to a much higher density (calculated as ~{densities['C']:.2f} g/cm³).\n"
                  f" - d) A 0.5-Earth-mass planet would be less compressed, leading to a lower density "
                  f"(calculated as ~{densities['D']:.2f} g/cm³).\n"
                  f"Therefore, planet (c) has the highest density.")
        return reason

# Run the check
print(check_exoplanet_density_answer())