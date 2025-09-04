import math

def check_exoplanet_density():
    """
    Checks which of the given exoplanet options has the highest density.

    The core physical principle is that for a planet of a given composition (like Earth),
    density increases with mass due to gravitational compression. We can model this with
    a mass-radius relationship R ∝ M^β, where β is an exponent.
    - For a constant density object, R ∝ M^(1/3), so β = 1/3.
    - For rocky planets, compression makes β < 1/3 (typically ~0.27).
    This implies that density ρ = M/V ∝ M/R^3 ∝ M/M^(3β) = M^(1-3β).
    Since β < 1/3, the exponent (1-3β) is positive, meaning density increases with mass.
    """

    # --- Define baseline and model parameters ---
    
    # Earth's average density in g/cm^3
    earth_density = 5.51
    
    # Exponent for the mass-density relationship: ρ_planet = ρ_earth * (M_planet/M_earth)^k
    # where k = 1 - 3*β. Using a typical β ≈ 0.27 for rocky planets, k ≈ 0.19.
    # Any small positive value for k will yield the same qualitative result.
    k = 0.19

    # --- Calculate density for each option ---

    # a) An Earth-mass and Earth-radius planet. This is our baseline.
    density_a = earth_density

    # b) A planet with 2 Earth masses and a density of approximately 5.5 g/cm^3.
    # The density is given in the problem statement.
    density_b = 5.5

    # c) A planet with the same composition as Earth but 5 times more massive.
    mass_ratio_c = 5.0
    density_c = earth_density * (mass_ratio_c ** k)

    # d) A planet with the same composition as Earth but half the mass of Earth.
    mass_ratio_d = 0.5
    density_d = earth_density * (mass_ratio_d ** k)

    # --- Verification ---

    # Store calculated densities in a dictionary for easy comparison.
    densities = {
        'a': density_a,
        'b': density_b,
        'c': density_c,
        'd': density_d
    }

    # The provided answer is 'C', which corresponds to option 'c'.
    llm_answer_option = 'c'

    # Find the option with the highest calculated density.
    highest_density_option = max(densities, key=densities.get)

    # Check if the LLM's answer matches the calculated result.
    if highest_density_option == llm_answer_option:
        return "Correct"
    else:
        # If the answer is wrong, provide a detailed reason.
        reason = (
            f"The provided answer '{llm_answer_option.upper()}' is incorrect.\n"
            f"The model calculates that option '{highest_density_option.upper()}' has the highest density.\n"
            "Here are the calculated densities based on physical principles:\n"
            f"  - a) Earth-like planet: {densities['a']:.2f} g/cm^3\n"
            f"  - b) 2 M_earth planet: {densities['b']:.2f} g/cm^3 (as given)\n"
            f"  - c) 5 M_earth planet: {densities['c']:.2f} g/cm^3 (denser due to compression)\n"
            f"  - d) 0.5 M_earth planet: {densities['d']:.2f} g/cm^3 (less dense due to lower compression)\n"
            f"The highest value is {densities[highest_density_option]:.2f} g/cm^3, which corresponds to option '{highest_density_option.upper()}'."
        )
        return reason

# Run the check
result = check_exoplanet_density()
print(result)