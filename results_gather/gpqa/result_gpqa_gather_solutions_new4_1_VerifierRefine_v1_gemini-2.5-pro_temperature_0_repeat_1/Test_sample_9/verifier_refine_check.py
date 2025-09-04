import math

def check_exoplanet_density_answer():
    """
    Checks the correctness of the LLM's answer about exoplanet density.

    The core principle is that for planets of the same composition, density increases
    with mass due to gravitational compression. This can be modeled with a mass-radius
    relationship R ∝ M^k, which leads to a density relationship ρ ∝ M^(1-3k).
    We use a typical exponent k=0.27 for rocky planets.
    """
    
    # The final answer provided by the LLM.
    llm_answer = "C"

    # --- Model Setup ---
    # Baseline density of Earth (g/cm^3)
    rho_earth = 5.51

    # Option a: Earth-mass and Earth-radius planet.
    density_a = rho_earth

    # Option b: 2 Earth masses and a density of approximately 5.5 g/cm^3.
    density_b = 5.5

    # For options c and d, we estimate density based on mass.
    # The exponent in the mass-density relationship: ρ ∝ M^exponent
    # Using k=0.27 for the mass-radius relation R ∝ M^k, the density exponent is 1 - 3*k
    mass_density_exponent = 1 - 3 * 0.27  # This is 0.19

    # Option c: 5 times more massive than Earth.
    mass_ratio_c = 5.0
    density_c = rho_earth * (mass_ratio_c ** mass_density_exponent)

    # Option d: Half the mass of Earth.
    mass_ratio_d = 0.5
    density_d = rho_earth * (mass_ratio_d ** mass_density_exponent)

    # --- Verification Logic ---
    
    # Store the descriptions and their calculated densities
    densities = {
        'a': density_a,
        'b': density_b,
        'c': density_c,
        'd': density_d
    }

    # Find the description corresponding to the highest density
    actual_densest_option = max(densities, key=densities.get)

    # Map the multiple-choice letters to the planet descriptions
    choice_to_description = {
        'A': 'a',
        'B': 'd',
        'C': 'c',
        'D': 'b'
    }
    
    # Get the description selected by the LLM
    llm_selected_description = choice_to_description.get(llm_answer)

    # Check if the LLM's choice matches the calculated result
    if llm_selected_description == actual_densest_option:
        return "Correct"
    else:
        # Format a detailed reason for the incorrect answer
        reason = f"The provided answer '{llm_answer}' is incorrect.\n"
        reason += "The reasoning should be that for planets of the same composition, higher mass leads to greater gravitational compression and thus higher density.\n"
        reason += "Based on this principle and quantitative estimates:\n"
        reason += f"  - Density(a) ≈ {densities['a']:.2f} g/cm³ (Earth baseline)\n"
        reason += f"  - Density(b) ≈ {densities['b']:.2f} g/cm³ (Given)\n"
        reason += f"  - Density(c) ≈ {densities['c']:.2f} g/cm³ (More massive, so denser)\n"
        reason += f"  - Density(d) ≈ {densities['d']:.2f} g/cm³ (Less massive, so less dense)\n"
        reason += f"The planet with the highest density is '{actual_densest_option}'. The correct multiple-choice option is the one that corresponds to description '{actual_densest_option}'.\n"
        reason += f"The LLM chose '{llm_answer}', which corresponds to description '{llm_selected_description}'."
        return reason

# Execute the check and print the result
print(check_exoplanet_density_answer())